"""
    plot_tds_flux(source; kwargs...) -> Figure

Plot desorption flux vs temperature for a TDS simulation.

`source` can be:
- A `String` path to a summary CSV written by `write_summary_csv`
- A `DataFrame` already loaded
- A `SimulationResult` (temperature must be in `result.summaries[:extra_timeseries]`
  under the key `"temperature_K"`)

# Keyword arguments
- `flux_vars` — variable names to plot fluxes for (default: all variables with flux columns).
- `temperature_column` — name of the temperature column in the CSV/DataFrame (default: `"temperature_K"`).
- `surface` — `:right`, `:left`, or `:both` (default: `:right`).
- `title` — plot title string (default: `"TDS Desorption Spectrum"`).
- `kwargs...` — forwarded to the Makie `lines!` call.

# Example
```julia
using Flopsy, CairoMakie

# From a summary CSV
fig = plot_tds_flux("tds_summary.csv")
save("tds_spectrum.png", fig)

# From a SimulationResult (temperature must be stored as extra_timeseries)
fig = plot_tds_flux(result)
```

Requires a Makie backend (`using CairoMakie` or `using GLMakie`).
"""
function plot_tds_flux(source; kwargs...)
    _makie_error("plot_tds_flux")
end


"""
    plot_spatial_snapshot(source; kwargs...) -> Figure

Plot the spatial profile of one or more variables at a single time point.

`source` can be:
- A `SimulationResult`
- A `String` path to an HDF5 file written by `write_field_output_hdf5`

# Keyword arguments
- `vars` — variable name(s) to plot (Symbol or Vector{Symbol}; default: all variables).
- `time_index` — saved-time index to plot (integer, or `:last`; default: `:last`).
- `group_fn` — optional grouping function `(var_names, data_matrix) -> Dict{String, Vector}`
  that aggregates variables (e.g. sum by trap type) before plotting.  Useful for Palioxis
  models with many occupancy-level DOFs.
- `title` — plot title (default: auto-generated from time).
- `xlabel` — x-axis label (default: `"Position (m)"`).
- `ylabel` — y-axis label (default: `"Concentration"`).

# Palioxis grouping helpers

For Palioxis-backed models with many trap occupancy levels, use a grouping function to
aggregate by trap type or by occupancy level before plotting:

```julia
# Sum all trap_H_* variables into a single "Total Trapped" series
group_fn = (names, data) -> begin
    mobile = Dict{String,Vector}(n => data[i, :] for (i, n) in enumerate(names) if startswith(n, "mobile"))
    trapped = sum(data[i, :] for (i, n) in enumerate(names) if startswith(n, "trap"); init=zeros(size(data, 2)))
    merge(mobile, Dict("Total Trapped" => trapped))
end

fig = plot_spatial_snapshot(result; group_fn=group_fn)
```

Requires a Makie backend.
"""
function plot_spatial_snapshot(source; kwargs...)
    _makie_error("plot_spatial_snapshot")
end


"""
    plot_spatial_evolution(source; kwargs...) -> Figure

Plot spatial profiles at multiple time points overlaid on a single axis, coloured by time.

`source` can be:
- A `SimulationResult`
- A `String` path to an HDF5 file

# Keyword arguments
- `vars` — variable name(s) to include (Symbol or Vector{Symbol}; default: all).
- `times` — time values to plot (matched to nearest saved time; default: 10 evenly-spaced
  snapshots across the full time range).
- `group_fn` — optional aggregation function (same interface as `plot_spatial_snapshot`).
- `colormap` — Makie colormap name (default: `"viridis"`).
- `title` — plot title.

Requires a Makie backend.
"""
function plot_spatial_evolution(source; kwargs...)
    _makie_error("plot_spatial_evolution")
end


"""
    record_spatial_video(source, output_path; kwargs...) -> output_path

Record an animated video of the spatial distribution evolving over time.

`source` can be:
- A `SimulationResult`
- A `String` path to an HDF5 file

# Keyword arguments
- `vars` — variable name(s) to animate (Symbol or Vector{Symbol}; default: all).
- `fps` — frames per second (default: `30`).
- `group_fn` — optional aggregation function.
- `title_fn` — function `t -> String` that generates a title from the current time
  (default: `t -> @sprintf("t = %.1f s", t)`).
- `ylims` — fixed y-axis limits `(ymin, ymax)`, or `nothing` for auto-scaling (default: `nothing`).

# Example
```julia
using Flopsy, CairoMakie

record_spatial_video(result, "spatial_evolution.mp4"; fps=24, vars=[:mobile_H])
```

Requires a Makie backend.
"""
function record_spatial_video(source, output_path; kwargs...)
    _makie_error("record_spatial_video")
end


# Internal helper -----------------------------------------------------------------

function _makie_error(fname::String)
    error("""
    $fname requires a Makie backend.  Load one before calling this function:

        using CairoMakie   # for file output (PNG, SVG, PDF, MP4)
        using GLMakie      # for interactive windows

    Then call $fname(...) again.
    """)
end
