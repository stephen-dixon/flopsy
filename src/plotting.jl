"""
    plot_tds_flux(source; kwargs...) -> Figure

Plot desorption flux vs temperature for a TDS simulation.

`source` can be a `String` path to a summary CSV, a `DataFrame`, or a `SimulationResult`.

# Keyword arguments
- `flux_vars` — variable names to plot (default: all flux columns).
- `temperature_column` — name of the temperature column (default: `"temperature_K"`).
- `surface` — `:right`, `:left`, or `:both` (default: `:right`).
- `title` — plot title (default: `"TDS Desorption Spectrum"`).

Requires a Makie backend (`using CairoMakie` or `using GLMakie`).
"""
function plot_tds_flux(source; kwargs...)
    _makie_error("plot_tds_flux")
end


"""
    plot_spatial_snapshot(source; kwargs...) -> Figure

Plot the spatial profile of one or more variables at a single time point.

`source` can be a `SimulationResult` or a `String` path to an HDF5 file.

# Keyword arguments
- `vars` — variable name(s) to plot (Symbol or Vector{Symbol}; default: all).
- `time_index` — saved-time index to plot (integer or `:last`; default: `:last`).
- `group_fn` — optional aggregation function `(names, data_matrix) -> Dict{String,Vector}`.
- `all_on_one` — plot all variables/groups on a single axis (default: `true`).
- `xscale` — `:log10` or `:identity` (default: `:log10`).
- `yscale` — `:log10` or `:identity` (default: `:log10`).
- `xmin` — clamp x values to this minimum before log scaling (default: `1e-20`).
- `ymin` — clamp y values to this minimum before log scaling (default: `1e-20`).
- `title` — plot title (default: auto from time).
- `xlabel` / `ylabel` — axis labels.

Requires a Makie backend.
"""
function plot_spatial_snapshot(source; kwargs...)
    _makie_error("plot_spatial_snapshot")
end


"""
    plot_spatial_evolution(source; kwargs...) -> Figure

Plot spatial profiles at multiple time points, coloured by time.

`source` can be a `SimulationResult` or a `String` path to an HDF5 file.

# Keyword arguments
- `vars` — variable name(s) (default: all).
- `times` — time values to plot (default: 10 evenly-spaced snapshots).
- `group_fn` — optional aggregation function.
- `all_on_one` — all variables on a single axis (default: `true`).
- `xscale` — `:log10` or `:identity` (default: `:log10`).
- `yscale` — `:log10` or `:identity` (default: `:log10`).
- `xmin` — minimum x clamp for log scaling (default: `1e-20`).
- `ymin` — minimum y clamp for log scaling (default: `1e-20`).
- `colormap` — Makie colormap (default: `:viridis`).
- `title` — plot title.

Requires a Makie backend.
"""
function plot_spatial_evolution(source; kwargs...)
    _makie_error("plot_spatial_evolution")
end


"""
    record_spatial_video(source, output_path; kwargs...) -> output_path

Record an animated video (e.g. MP4) of the spatial distribution per variable
evolving over time.  Variables are plotted in separate panels.

`source` can be a `SimulationResult` or a `String` path to an HDF5 file.

# Keyword arguments
- `vars` — variable name(s) to animate (default: all).
- `fps` — frames per second (default: `30`).
- `group_fn` — optional aggregation function.
- `title_fn` — `t -> String` (default: `t -> "t = %.3g s"`).
- `ylims` — fixed y limits or `nothing` for auto-scaling.

Requires a Makie backend.
"""
function record_spatial_video(source, output_path; kwargs...)
    _makie_error("record_spatial_video")
end


"""
    record_spatial_animation(source, output_path; kwargs...) -> output_path

Record an animated video of the spatial distribution with all variables/groups
on a single log-log axis, evolving over every saved time step.

`source` can be a `SimulationResult` or a `String` path to an HDF5 file.

# Keyword arguments
- `vars` — variable name(s) (default: all).
- `fps` — frames per second (default: `30`).
- `group_fn` — optional aggregation function.
- `title_fn` — `t -> String` title generator (default: `t -> "t = %.3g s"`).
- `xscale` — `:log10` or `:identity` (default: `:log10`).
- `yscale` — `:log10` or `:identity` (default: `:log10`).
- `xmin` — minimum x clamp for log scaling (default: `1e-20`).
- `ymin` — minimum y clamp for log scaling (default: `1e-20`).
- `ylims` — fixed y limits after clamping, or `nothing` for auto.
- `colormap` — colormap for line colours cycling over variables (default: `:tab10`).

Requires a Makie backend (`using CairoMakie` for MP4 output).

# Example
```julia
using Flopsy, CairoMakie
record_spatial_animation(result, "spatial.mp4"; fps=24, group_fn=my_group_fn)
```
"""
function record_spatial_animation(source, output_path; kwargs...)
    _makie_error("record_spatial_animation")
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
