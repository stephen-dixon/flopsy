__precompile__(false)

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
            fields[name] = read(h5["fields"][name])
        end
    end
    return t, x, dx, fields
end

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

function _apply_grouping(fields::Dict{String,<:AbstractMatrix}, it::Int, group_fn)
    names = collect(keys(fields))
    nx = size(first(values(fields)), 2)
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

function _filter_vars(fields::Dict{String,<:AbstractMatrix}, vars)
    vars === nothing && return fields
    wanted = vars isa Symbol ? [string(vars)] : string.(vars)
    return Dict(k => v for (k, v) in fields if k ∈ wanted)
end

function _nearest_time_index(t::AbstractVector, t_req)
    t_req === :last && return length(t)
    return argmin(abs.(t .- Float64(t_req)))
end

# Clamp values and apply log scaling to x or y vectors.
_clamp_log(v::AbstractVector, min_val::Real, scale::Symbol) =
    scale === :log10 ? max.(v, min_val) : v

_scale_sym(s::Symbol) = s === :log10 ? Makie.log10 : Makie.identity


# ---------------------------------------------------------------------------
# plot_tds_flux
# ---------------------------------------------------------------------------

function Flopsy.plot_tds_flux(source; kwargs...)
    _plot_tds_flux_impl(source; kwargs...)
end

function _plot_tds_flux_impl(csv_path::AbstractString;
    flux_vars          = nothing,
    temperature_column = "temperature_K",
    surface            = :right,
    title              = "TDS Desorption Spectrum",
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

    if T_col ∈ all_cols
        T_vec = Float64.(df[!, T_col])
    elseif "time" ∈ all_cols
        @warn "Temperature column '$T_col' not found; using time as x-axis"
        T_vec = Float64.(df[!, "time"])
        T_col = "time"
    else
        error("Cannot find a suitable x-axis column in the DataFrame.")
    end

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
    all_on_one = true,
    xscale     = :log10,
    yscale     = :log10,
    xmin       = 1e-20,
    ymin       = 1e-20,
    title      = nothing,
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _extract_result_fields(result)
    return _spatial_snapshot_plot(x, t, fields, vars, time_index, group_fn,
                                   all_on_one, xscale, yscale, xmin, ymin,
                                   title, xlabel, ylabel; kwargs...)
end

function _plot_spatial_snapshot_impl(hdf5_path::AbstractString;
    vars       = nothing,
    time_index = :last,
    group_fn   = nothing,
    all_on_one = true,
    xscale     = :log10,
    yscale     = :log10,
    xmin       = 1e-20,
    ymin       = 1e-20,
    title      = nothing,
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _load_hdf5_fields(hdf5_path)
    return _spatial_snapshot_plot(x, t, fields, vars, time_index, group_fn,
                                   all_on_one, xscale, yscale, xmin, ymin,
                                   title, xlabel, ylabel; kwargs...)
end

function _spatial_snapshot_plot(x, t, fields, vars, time_index, group_fn,
                                  all_on_one, xscale, yscale, xmin, ymin,
                                  title_in, xlabel, ylabel; kwargs...)
    fields = _filter_vars(fields, vars)
    it     = _nearest_time_index(t, time_index)
    t_val  = t[it]

    title  = title_in !== nothing ? title_in : @sprintf("t = %.3g s", t_val)
    grouped = _apply_grouping(fields, it, group_fn)

    x_plot = _clamp_log(x, xmin, xscale)

    if all_on_one || length(grouped) == 1
        fig = Figure(size = (720, 440))
        ax  = Axis(fig[1, 1];
                   xlabel    = xlabel,
                   ylabel    = ylabel,
                   title     = title,
                   xscale    = _scale_sym(xscale),
                   yscale    = _scale_sym(yscale))

        for (name, profile) in sort(collect(grouped); by=first)
            y_plot = _clamp_log(profile, ymin, yscale)
            lines!(ax, x_plot, y_plot; label = name, kwargs...)
        end

        length(grouped) > 1 && axislegend(ax; position = :rt)
    else
        var_names = collect(sort(collect(keys(grouped))))
        n_vars    = length(var_names)
        fig = Figure(size = (720, 320 * n_vars))

        for (ivar, vname) in enumerate(var_names)
            ax = Axis(fig[ivar, 1];
                      xlabel = ivar == n_vars ? xlabel : "",
                      ylabel = ylabel,
                      title  = ivar == 1 ? "$title — $vname" : vname,
                      xscale = _scale_sym(xscale),
                      yscale = _scale_sym(yscale))
            y_plot = _clamp_log(grouped[vname], ymin, yscale)
            lines!(ax, x_plot, y_plot; kwargs...)
        end
    end

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
    all_on_one = true,
    xscale     = :log10,
    yscale     = :log10,
    xmin       = 1e-20,
    ymin       = 1e-20,
    colormap   = :viridis,
    title      = "Spatial Evolution",
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _extract_result_fields(result)
    return _spatial_evolution_plot(x, t, fields, vars, times, group_fn,
                                    all_on_one, xscale, yscale, xmin, ymin,
                                    colormap, title, xlabel, ylabel; kwargs...)
end

function _plot_spatial_evolution_impl(hdf5_path::AbstractString;
    vars       = nothing,
    times      = nothing,
    group_fn   = nothing,
    all_on_one = true,
    xscale     = :log10,
    yscale     = :log10,
    xmin       = 1e-20,
    ymin       = 1e-20,
    colormap   = :viridis,
    title      = "Spatial Evolution",
    xlabel     = "Position (m)",
    ylabel     = "Concentration",
    kwargs...
)
    t, x, dx, fields = _load_hdf5_fields(hdf5_path)
    return _spatial_evolution_plot(x, t, fields, vars, times, group_fn,
                                    all_on_one, xscale, yscale, xmin, ymin,
                                    colormap, title, xlabel, ylabel; kwargs...)
end

function _spatial_evolution_plot(x, t, fields, vars, times_req, group_fn,
                                   all_on_one, xscale, yscale, xmin, ymin,
                                   colormap, title, xlabel, ylabel; kwargs...)
    fields = _filter_vars(fields, vars)
    nt     = length(t)

    if times_req === nothing
        n_snap   = min(10, nt)
        snap_idx = round.(Int, LinRange(1, nt, n_snap))
    else
        snap_idx = [_nearest_time_index(t, tv) for tv in times_req]
    end
    snap_idx = unique(snap_idx)

    cmap   = to_colormap(colormap)
    colors = [cmap[round(Int, 1 + (i-1)/(max(length(snap_idx)-1,1)) * (length(cmap)-1))]
              for i in 1:length(snap_idx)]

    x_plot = _clamp_log(x, xmin, xscale)

    # Determine variable/group names at first snap to know what we're plotting.
    sample_grouped = _apply_grouping(fields, snap_idx[1], group_fn)
    var_names = sort(collect(keys(sample_grouped)))
    n_vars    = length(var_names)

    if all_on_one || n_vars == 1
        fig = Figure(size = (720, 440))
        ax  = Axis(fig[1, 1];
                   xlabel  = xlabel,
                   ylabel  = ylabel,
                   title   = title,
                   xscale  = _scale_sym(xscale),
                   yscale  = _scale_sym(yscale))

        for (k, it) in enumerate(snap_idx)
            grouped = _apply_grouping(fields, it, group_fn)
            for vname in var_names
                profile = haskey(grouped, vname) ? grouped[vname] : fields[vname][it, :]
                y_plot  = _clamp_log(profile, ymin, yscale)
                lines!(ax, x_plot, y_plot;
                       color = colors[k],
                       label = @sprintf("t=%.2g", t[it]),
                       kwargs...)
            end
        end

        Colorbar(fig[1, end+1];
                 colormap = colormap,
                 limits   = (t[snap_idx[1]], t[snap_idx[end]]),
                 label    = "Time (s)")
    else
        fig = Figure(size = (720, 320 * n_vars))

        for (ivar, vname) in enumerate(var_names)
            ax = Axis(fig[ivar, 1];
                      xlabel = ivar == n_vars ? xlabel : "",
                      ylabel = ylabel,
                      title  = ivar == 1 ? "$title — $vname" : vname,
                      xscale = _scale_sym(xscale),
                      yscale = _scale_sym(yscale))

            for (k, it) in enumerate(snap_idx)
                grouped = _apply_grouping(fields, it, group_fn)
                profile = haskey(grouped, vname) ? grouped[vname] : fields[vname][it, :]
                y_plot  = _clamp_log(profile, ymin, yscale)
                lines!(ax, x_plot, y_plot;
                       color = colors[k],
                       label = @sprintf("t=%.2g", t[it]),
                       kwargs...)
            end
        end

        Colorbar(fig[:, end+1];
                 colormap = colormap,
                 limits   = (t[snap_idx[1]], t[snap_idx[end]]),
                 label    = "Time (s)")
    end

    return fig
end


# ---------------------------------------------------------------------------
# record_spatial_video  (per-variable panels, linear scale by default)
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

    if ylims === nothing
        ymin_val = minimum(minimum(fields[v]) for v in var_names)
        ymax_val = maximum(maximum(fields[v]) for v in var_names)
        pad      = 0.05 * (ymax_val - ymin_val + 1e-30)
        y_limits = (ymin_val - pad, ymax_val + pad)
    else
        y_limits = ylims
    end

    fig    = Figure(size = (720, 320 * max(1, n_vars)))
    axes   = [Axis(fig[i, 1];
                   xlabel = i == n_vars ? "Position (m)" : "",
                   ylabel = "Concentration",
                   title  = var_names[i],
                   limits = (nothing, y_limits))
              for i in 1:n_vars]

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


# ---------------------------------------------------------------------------
# record_spatial_animation  (all variables on one log-log axis)
# ---------------------------------------------------------------------------

function Flopsy.record_spatial_animation(source, output_path; kwargs...)
    _record_spatial_animation_impl(source, output_path; kwargs...)
end

function _record_spatial_animation_impl(result::SimulationResult, output_path;
    vars      = nothing,
    fps       = 30,
    group_fn  = nothing,
    title_fn  = nothing,
    xscale    = :log10,
    yscale    = :log10,
    xmin      = 1e-20,
    ymin      = 1e-20,
    ylims     = nothing,
    colormap  = :tab10,
    kwargs...
)
    t, x, dx, fields = _extract_result_fields(result)
    return _do_animate(x, t, fields, output_path, vars, fps, group_fn, title_fn,
                        xscale, yscale, xmin, ymin, ylims, colormap; kwargs...)
end

function _record_spatial_animation_impl(hdf5_path::AbstractString, output_path;
    vars      = nothing,
    fps       = 30,
    group_fn  = nothing,
    title_fn  = nothing,
    xscale    = :log10,
    yscale    = :log10,
    xmin      = 1e-20,
    ymin      = 1e-20,
    ylims     = nothing,
    colormap  = :tab10,
    kwargs...
)
    t, x, dx, fields = _load_hdf5_fields(hdf5_path)
    return _do_animate(x, t, fields, output_path, vars, fps, group_fn, title_fn,
                        xscale, yscale, xmin, ymin, ylims, colormap; kwargs...)
end

function _do_animate(x, t, fields, output_path, vars, fps, group_fn, title_fn,
                      xscale, yscale, xmin, ymin, ylims, colormap; kwargs...)
    fields = _filter_vars(fields, vars)
    nt     = length(t)

    default_title_fn = tv -> @sprintf("t = %.3g s", tv)
    tfn = title_fn !== nothing ? title_fn : default_title_fn

    # Determine group names from first time step.
    sample = _apply_grouping(fields, 1, group_fn)
    var_names = sort(collect(keys(sample)))
    n_vars    = length(var_names)

    cmap   = to_colormap(colormap)
    colors = [cmap[round(Int, 1 + (i-1)/(max(n_vars-1,1)) * (length(cmap)-1))]
              for i in 1:n_vars]

    x_plot = _clamp_log(x, xmin, xscale)

    # Determine y limits across all times and variables.
    if ylims === nothing
        all_vals = Float64[]
        for it in 1:nt
            grouped = _apply_grouping(fields, it, group_fn)
            for vname in var_names
                profile = haskey(grouped, vname) ? grouped[vname] : fields[vname][it, :]
                append!(all_vals, _clamp_log(profile, ymin, yscale))
            end
        end
        ymin_v = minimum(all_vals)
        ymax_v = maximum(all_vals)
        pad    = 0.0  # log scale — don't add linear padding
        y_limits = (ymin_v, ymax_v)
    else
        y_limits = ylims
    end

    fig = Figure(size = (760, 480))
    ax  = Axis(fig[1, 1];
               xlabel = "Position (m)",
               ylabel = "Concentration",
               xscale = _scale_sym(xscale),
               yscale = _scale_sym(yscale),
               limits = (nothing, y_limits))

    frame_idx = Observable(1)
    title_obs = @lift tfn(t[$frame_idx])
    Label(fig[0, 1], title_obs; fontsize = 14)

    for (i, vname) in enumerate(var_names)
        profile_obs = @lift begin
            grouped = _apply_grouping(fields, $frame_idx, group_fn)
            raw = haskey(grouped, vname) ? Float64.(grouped[vname]) : Float64.(fields[vname][$frame_idx, :])
            _clamp_log(raw, ymin, yscale)
        end
        lines!(ax, x_plot, profile_obs;
               color = colors[i],
               label = vname,
               kwargs...)
    end

    n_vars > 1 && axislegend(ax; position = :rt)

    Makie.record(fig, output_path, 1:nt; framerate = fps) do it
        frame_idx[] = it
    end

    return output_path
end

end # module MakieExt
