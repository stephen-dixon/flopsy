module MakieExt

using Flopsy
using Makie
using CSV
using DataFrames
using HDF5
import Printf: @sprintf

# ---------------------------------------------------------------------------
# Internal data-loading helpers
# ---------------------------------------------------------------------------

"""Load (t, x, dx, fields::Dict{String,Matrix}) from an HDF5 path."""
function _load_hdf5_fields(path::AbstractString)
    t = Float64[]
    x = Float64[]
    dx = 0.0
    fields = Dict{String,Matrix{Float64}}()
    h5open(path, "r") do h5
        t  = read(h5["time"])
        x  = read(h5["mesh/x"])
        dx = read(h5["mesh/dx"])
        for name in keys(h5["fields"])
            fields[name] = read(h5["fields"][name])   # (nt, nx)
        end
    end
    return t, x, dx, fields
end

"""Extract (t, x, dx, fields::Dict{String,Matrix}) from a SimulationResult."""
function _extract_result_fields(result::SimulationResult)
    model  = result.model
    layout = model.layout
    nx     = model.context.nx
    x      = model.context.mesh.x
    dx     = model.context.mesh.dx
    t      = result.solution.t

    fields = Dict{String,Matrix{Float64}}()
    for name in variable_names(layout)
        fields[string(name)] = variable_timeseries(result, name)
    end
    return collect(Float64, t), x, dx, fields
end

"""Apply optional group_fn to (var_names, fields) → Dict{String, Vector} at time index it."""
function _apply_grouping(fields::Dict{String,<:AbstractMatrix}, it::Int, group_fn)
    names = collect(keys(fields))
    nt, nx = size(first(values(fields)))
    # Build a matrix view at this time index
    data = Matrix{Float64}(undef, length(names), nx)
    for (i, n) in enumerate(names)
        data[i, :] .= fields[n][it, :]
    end
    if group_fn === nothing
        return Dict(n => fields[n][it, :] for n in names)
    else
        return group_fn(names, data)
    end
end

"""Select a subset of variables from a fields dict."""
function _filter_vars(fields::Dict{String,<:AbstractMatrix}, vars)
    vars === nothing && return fields
    wanted = vars isa Symbol ? [string(vars)] : string.(vars)
    return Dict(k => v for (k, v) in fields if k ∈ wanted)
end

"""Find the saved-time index nearest to a requested time value."""
function _nearest_time_index(t::AbstractVector, t_req)
    t_req === :last && return length(t)
    return argmin(abs.(t .- Float64(t_req)))
end

# ---------------------------------------------------------------------------
# plot_tds_flux
# ---------------------------------------------------------------------------

function Flopsy.plot_tds_flux(source; kwargs...)
    _plot_tds_flux_impl(source; kwargs...)
end

function _plot_tds_flux_impl(csv_path::AbstractString;
    flux_vars           = nothing,
    temperature_column  = "temperature_K",
    surface             = :right,
    title               = "TDS Desorption Spectrum",
    kwargs...
)
    df = CSV.read(csv_path, DataFrame)
    return _tds_flux_from_df(df, flux_vars, temperature_column, surface, title; kwargs...)
end

function _plot_tds_flux_impl(df::DataFrame;
    flux_vars          = nothing,
    temperature_column = "temperature_K",
    surface            = :right,
    title              = "TDS Desorption Spectrum",
    kwargs...
)
    return _tds_flux_from_df(df, flux_vars, temperature_column, surface, title; kwargs...)
end

function _plot_tds_flux_impl(result::SimulationResult;
    flux_vars          = nothing,
    temperature_column = "temperature_K",
    surface            = :right,
    title              = "TDS Desorption Spectrum",
    kwargs...
)
    df = build_summary_dataframe(result)
    return _tds_flux_from_df(df, flux_vars, temperature_column, surface, title; kwargs...)
end

function _tds_flux_from_df(df, flux_vars, T_col, surface, title; kwargs...)
    all_cols = string.(names(df))

    # Temperature axis
    if T_col ∈ all_cols
        T_vec = Float64.(df[!, T_col])
    elseif "time" ∈ all_cols
        @warn "Temperature column '$T_col' not found; using time as x-axis"
        T_vec = Float64.(df[!, "time"])
        T_col = "time"
    else
        error("Cannot find a suitable x-axis column in the DataFrame.")
    end

    # Determine flux columns to plot
    prefix = surface === :left  ? "left_flux_" :
             surface === :right ? "right_flux_" : nothing

    if prefix !== nothing
        flux_cols = filter(c -> startswith(c, prefix), all_cols)
    else
        flux_cols = filter(c -> startswith(c, "left_flux_") || startswith(c, "right_flux_"), all_cols)
    end

    if flux_vars !== nothing
        wanted = Set(string.(flux_vars isa Symbol ? [flux_vars] : flux_vars))
        flux_cols = filter(c -> any(endswith(c, w) for w in wanted), flux_cols)
    end

    isempty(flux_cols) && error("No flux columns found (surface=:$surface, flux_vars=$flux_vars)")

    fig = Figure(size = (800, 450))
    ax  = Axis(fig[1, 1];
               xlabel = T_col == "time" ? "Time (s)" : "Temperature (K)",
               ylabel = "Desorption flux",
               title  = title)

    for col in flux_cols
        y = Float64.(df[!, col])
        label = replace(col, "right_flux_" => "", "left_flux_" => "")
        lines!(ax, T_vec, y; label = label, kwargs...)
    end

    length(flux_cols) > 1 && axislegend(ax; position = :rt)

    return fig
end

# ---------------------------------------------------------------------------
# plot_spatial_snapshot
# ---------------------------------------------------------------------------

function Flopsy.plot_spatial_snapshot(source; kwargs...)
    _plot_spatial_snapshot_impl(source; kwargs...)
end

function _plot_spatial_snapshot_impl(result::SimulationResult;
    vars       = nothing,
    time_index = :last,
    group_fn   = nothing,
    title      = nothing,
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _extract_result_fields(result)
    return _spatial_snapshot_plot(x, t, fields, vars, time_index, group_fn, title, xlabel, ylabel; kwargs...)
end

function _plot_spatial_snapshot_impl(hdf5_path::AbstractString;
    vars       = nothing,
    time_index = :last,
    group_fn   = nothing,
    title      = nothing,
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _load_hdf5_fields(hdf5_path)
    return _spatial_snapshot_plot(x, t, fields, vars, time_index, group_fn, title, xlabel, ylabel; kwargs...)
end

function _spatial_snapshot_plot(x, t, fields, vars, time_index, group_fn, title_in, xlabel, ylabel; kwargs...)
    fields = _filter_vars(fields, vars)
    it     = _nearest_time_index(t, time_index)
    t_val  = t[it]

    title  = title_in !== nothing ? title_in : @sprintf("t = %.3g s", t_val)

    fig = Figure(size = (700, 420))
    ax  = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, title = title)

    grouped = _apply_grouping(fields, it, group_fn)
    for (name, profile) in grouped
        lines!(ax, x, profile; label = name, kwargs...)
    end

    length(grouped) > 1 && axislegend(ax; position = :rt)

    return fig
end

# ---------------------------------------------------------------------------
# plot_spatial_evolution
# ---------------------------------------------------------------------------

function Flopsy.plot_spatial_evolution(source; kwargs...)
    _plot_spatial_evolution_impl(source; kwargs...)
end

function _plot_spatial_evolution_impl(result::SimulationResult;
    vars       = nothing,
    times      = nothing,
    group_fn   = nothing,
    colormap   = :viridis,
    title      = "Spatial Evolution",
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _extract_result_fields(result)
    return _spatial_evolution_plot(x, t, fields, vars, times, group_fn, colormap, title, xlabel, ylabel; kwargs...)
end

function _plot_spatial_evolution_impl(hdf5_path::AbstractString;
    vars       = nothing,
    times      = nothing,
    group_fn   = nothing,
    colormap   = :viridis,
    title      = "Spatial Evolution",
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _load_hdf5_fields(hdf5_path)
    return _spatial_evolution_plot(x, t, fields, vars, times, group_fn, colormap, title, xlabel, ylabel; kwargs...)
end

function _spatial_evolution_plot(x, t, fields, vars, times_req, group_fn, colormap, title, xlabel, ylabel; kwargs...)
    fields = _filter_vars(fields, vars)
    nt     = length(t)

    # Choose time indices to plot
    if times_req === nothing
        n_snap = min(10, nt)
        snap_idx = round.(Int, LinRange(1, nt, n_snap))
    else
        snap_idx = [_nearest_time_index(t, tv) for tv in times_req]
    end
    snap_idx = unique(snap_idx)

    var_names = collect(keys(fields))
    n_vars    = length(var_names)

    fig = Figure(size = (700, 420 * max(1, n_vars)))

    cmap   = to_colormap(colormap)
    colors = [cmap[round(Int, 1 + (i-1)/(max(length(snap_idx)-1, 1)) * (length(cmap)-1))]
              for i in 1:length(snap_idx)]

    for (ivar, vname) in enumerate(var_names)
        ax = Axis(fig[ivar, 1];
                  xlabel = ivar == n_vars ? xlabel : "",
                  ylabel = ylabel,
                  title  = ivar == 1 ? "$title — $vname" : vname)

        for (k, it) in enumerate(snap_idx)
            grouped = _apply_grouping(fields, it, group_fn)
            profile = haskey(grouped, vname) ? grouped[vname] : fields[vname][it, :]
            lines!(ax, x, profile; color = colors[k], label = @sprintf("t=%.2g", t[it]), kwargs...)
        end
    end

    # Shared colorbar to indicate time
    Colorbar(fig[:, end+1];
             colormap = colormap,
             limits   = (t[snap_idx[1]], t[snap_idx[end]]),
             label    = "Time (s)")

    return fig
end

# ---------------------------------------------------------------------------
# record_spatial_video
# ---------------------------------------------------------------------------

function Flopsy.record_spatial_video(source, output_path; kwargs...)
    _record_spatial_video_impl(source, output_path; kwargs...)
end

function _record_spatial_video_impl(result::SimulationResult, output_path;
    vars     = nothing,
    fps      = 30,
    group_fn = nothing,
    title_fn = nothing,
    ylims    = nothing,
    kwargs...
)
    t, x, dx, fields = _extract_result_fields(result)
    return _do_record(x, t, fields, output_path, vars, fps, group_fn, title_fn, ylims; kwargs...)
end

function _record_spatial_video_impl(hdf5_path::AbstractString, output_path;
    vars     = nothing,
    fps      = 30,
    group_fn = nothing,
    title_fn = nothing,
    ylims    = nothing,
    kwargs...
)
    t, x, dx, fields = _load_hdf5_fields(hdf5_path)
    return _do_record(x, t, fields, output_path, vars, fps, group_fn, title_fn, ylims; kwargs...)
end

function _do_record(x, t, fields, output_path, vars, fps, group_fn, title_fn, ylims; kwargs...)
    fields = _filter_vars(fields, vars)
    nt     = length(t)
    default_title_fn = tv -> @sprintf("t = %.3g s", tv)
    tfn    = title_fn !== nothing ? title_fn : default_title_fn

    var_names = collect(keys(fields))
    n_vars    = length(var_names)

    # Pre-compute y-axis limits across all times if not provided
    if ylims === nothing
        ymin = minimum(minimum(fields[v]) for v in var_names)
        ymax = maximum(maximum(fields[v]) for v in var_names)
        pad  = 0.05 * (ymax - ymin + 1e-30)
        y_limits = (ymin - pad, ymax + pad)
    else
        y_limits = ylims
    end

    fig    = Figure(size = (700, 420 * max(1, n_vars)))
    axes   = [Axis(fig[i, 1];
                   xlabel = i == n_vars ? "Position (m)" : "",
                   ylabel = "Concentration",
                   title  = var_names[i],
                   limits = (nothing, y_limits))
              for i in 1:n_vars]

    # Observable for time step
    frame_idx = Observable(1)
    title_obs = @lift tfn(t[$frame_idx])

    Label(fig[0, 1], title_obs; fontsize = 14)

    for (ivar, vname) in enumerate(var_names)
        profile_obs = @lift begin
            grouped = _apply_grouping(fields, $frame_idx, group_fn)
            haskey(grouped, vname) ? Float64.(grouped[vname]) : Float64.(fields[vname][$frame_idx, :])
        end
        lines!(axes[ivar], x, profile_obs; kwargs...)
    end

    Makie.record(fig, output_path, 1:nt; framerate = fps) do it
        frame_idx[] = it
    end

    return output_path
end

end # module MakieExt
