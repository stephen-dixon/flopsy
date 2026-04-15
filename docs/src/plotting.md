# Plotting

Flopsy provides optional plotting helpers via a [Makie](https://docs.makie.org/) package
extension.  The functions are exported from `Flopsy` but require a Makie backend to be
loaded before calling them.

## Setup

Install and load a Makie backend:

```julia
# Static file output (PNG, SVG, PDF, MP4) — recommended for scripts
using CairoMakie

# Interactive windows — recommended for interactive sessions
using GLMakie
```

Once a backend is loaded, the Flopsy plotting functions become available.

---

## TDS Desorption Spectrum

`plot_tds_flux` plots surface desorption flux vs temperature (or time) from a TDS summary.

### From a summary CSV

```julia
using Flopsy, CairoMakie

fig = plot_tds_flux("tds_summary.csv")
save("tds_spectrum.png", fig)
```

The CSV must have been written by `write_summary_csv`.  By default the function looks for
a `temperature_K` column; if absent it falls back to `time`.

### From a SimulationResult

```julia
fig = plot_tds_flux(result)
save("tds_spectrum.png", fig)
```

The `result` must have `"temperature_K"` in `result.summaries[:extra_timeseries]` (see the
TDS example in `examples/tds_palioxis.jl`).

### Options

```julia
fig = plot_tds_flux("tds_summary.csv";
    flux_vars          = [:mobile_H],      # restrict to specific variables
    surface            = :right,           # :left, :right, or :both
    temperature_column = "temperature_K",  # column name for x-axis
    title              = "TDS Spectrum",
)
```

---

## Spatial Snapshot

`plot_spatial_snapshot` plots the concentration profile at a single saved time.

### From a SimulationResult

```julia
using Flopsy, CairoMakie

fig = plot_spatial_snapshot(result)                       # final time
fig = plot_spatial_snapshot(result; time_index = 50)      # 50th saved time
save("snapshot.png", fig)
```

### From an HDF5 file

```julia
fig = plot_spatial_snapshot("tds_fields.h5"; time_index = :last)
```

### Options

```julia
fig = plot_spatial_snapshot(result;
    vars       = [:mobile_H, :trap_H_1],   # select subset of variables
    time_index = :last,                    # integer index or :last
    xlabel     = "Depth (m)",
    ylabel     = "Atom fraction",
    title      = "H profile at end of TDS",
)
```

### Grouping variables (Palioxis models)

Palioxis models with many occupancy-level DOFs can become cluttered.  Use `group_fn` to
aggregate variables before plotting.

```julia
# Sum all trap occupancy levels into a single "Total Trapped" series
group_fn = (names, data) -> begin
    mobile  = Dict(n => data[i, :] for (i, n) in enumerate(names) if startswith(n, "mobile"))
    trapped = sum(data[i, :] for (i, n) in enumerate(names)
                  if startswith(n, "trap"); init = zeros(size(data, 2)))
    merge(mobile, Dict("Total Trapped" => trapped))
end

fig = plot_spatial_snapshot(result; group_fn = group_fn)
```

The `group_fn` signature is:

```julia
group_fn(var_names::Vector{String}, data::Matrix{Float64}) -> Dict{String, Vector{Float64}}
```

where `data` has shape `(nvars, nx)` at the requested time index.

---

## Spatial Evolution

`plot_spatial_evolution` overlays profiles from multiple time snapshots on a single axis,
coloured by time.

### Basic usage

```julia
fig = plot_spatial_evolution(result)              # 10 evenly-spaced snapshots
fig = plot_spatial_evolution("tds_fields.h5")    # from HDF5

save("evolution.png", fig)
```

### Custom snapshot times

```julia
fig = plot_spatial_evolution(result;
    times    = [0.0, 100.0, 500.0, 1000.0, 3000.0],  # requested times (nearest saved used)
    vars     = [:mobile_H],
    colormap = :plasma,
    title    = "Mobile H evolution",
)
```

---

## Spatial Video

`record_spatial_video` renders an animation of the spatial distribution evolving over all
saved time steps and writes it to a video file.

```julia
using Flopsy, CairoMakie

# From a SimulationResult
record_spatial_video(result, "spatial_evolution.mp4"; fps = 24)

# From an HDF5 file
record_spatial_video("tds_fields.h5", "spatial_evolution.mp4"; fps = 30)
```

### Options

```julia
record_spatial_video(result, "out.mp4";
    vars     = [:mobile_H],             # restrict to specific variables
    fps      = 30,                      # frames per second
    ylims    = (0.0, 0.1),              # fixed y-axis (default: auto-scale across all times)
    title_fn = t -> "T = $(round(300 + t/60; digits=1)) K",   # custom title per frame
    group_fn = group_fn,                # optional variable aggregation
)
```

!!! tip
    For best results with MP4 output, use CairoMakie.  GLMakie records to `.mkv` or
    `.webm` by default; the output format is inferred from the file extension.

---

## API Reference

See the [API Reference](api.md#Plotting) for full docstring listings.
