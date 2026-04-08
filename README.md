# flopsy

![Logo](./docs/flopsy.png)

A Julia package for stiff 1D reaction-diffusion problems, targeting hydrogen transport
through materials.  Composable operators, multiple time-integration formulations, and
integration with the Palioxis trapping library.

---

## Quickstart

### Build a trapping model and solve it

```julia
using Flopsy
using OrdinaryDiffEq

# 1 mm domain, 60 nodes
mesh = Mesh1D(1e-3, 60)
nx   = length(mesh.x)

# Simple mobile + trap system
model = build_trapping_model(
    mesh                  = mesh,
    k_trap                = 5e-3,
    k_detrap              = 0.05,
    diffusion_coefficient = 1e-7,
)

# Gaussian initial mobile profile; traps empty
x   = mesh.x
u0  = zeros(2 * nx)
for ix in 1:nx
    u0[(ix-1)*2 + 1] = 0.05 * exp(-(x[ix] - x[end]/2)^2 / (2*(x[end]/8)^2))
end

tspan  = (0.0, 200.0)
saveat = collect(0.0:2.0:200.0)

config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm   = Rodas5P(),
    abstol      = 1e-9,
    reltol      = 1e-7,
    saveat      = saveat,
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
write_summary_csv(result, "summary.csv")
write_field_output_hdf5(result, "fields.h5")
```

---

## Switching Formulations

All formulations share the same model, initial condition, and output helpers.
Swap the `SolverConfig` to change the numerical technique:

```julia
using Flopsy, OrdinaryDiffEq, ADTypes

# ── Monolithic stiff ODE (default, uses analytic sparse Jacobian) ──────────
config_unsplit = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm   = Rodas5P(),
    abstol = 1e-9, reltol = 1e-7, saveat = saveat,
)

# ── IMEX: diffusion implicit, reaction explicit ────────────────────────────
config_imex = SolverConfig(
    formulation = IMEXFormulation(),
    algorithm   = KenCarp4(),
    abstol = 1e-9, reltol = 1e-7, saveat = saveat,
)

# ── Strang operator splitting (2nd-order) ─────────────────────────────────
config_strang = SolverConfig(
    formulation = SplitFormulation(StrangSplit()),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    dt          = 2.0,          # macro-step size (seconds)
    abstol = 1e-9, reltol = 1e-7, saveat = saveat,
)

# ── Lie operator splitting (1st-order) ────────────────────────────────────
config_lie = SolverConfig(
    formulation = SplitFormulation(LieSplit()),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    dt          = 2.0,
    abstol = 1e-9, reltol = 1e-7, saveat = saveat,
)

# ── DAE with quasi-static trapping constraint ─────────────────────────────
config_dae = SolverConfig(
    formulation = ResidualFormulation(),
    algorithm   = Rodas5P(),
    abstol = 1e-9, reltol = 1e-7, saveat = saveat,
)

# Run each with the same model and initial condition
results = map([config_unsplit, config_imex, config_strang, config_lie, config_dae]) do cfg
    wrap_result(model, solve_problem(model, u0, tspan, cfg), cfg)
end
```

See `examples/formulation_comparison.jl` for a full benchmark comparison.

---

## TDS Simulation with Palioxis

```julia
using Flopsy, Palioxis, OrdinaryDiffEq, ADTypes

Palioxis.init(ENV["PALIOXIS_ROOT"])
pal = MultipleDefectModel(ENV["PALIOXIS_XML"])

mesh    = Mesh1D(1e-3, 200)
defects = fill(1e-3, get_n_trap_types(pal), 200)

# Piecewise temperature: load at 400 K, rest at 300 K, ramp 300→1200 K at 10 K/min
T_load, T_rest, ramp_rate = 400.0, 300.0, 10.0/60.0
t_load, t_rest_end = 86400.0, 93600.0
T_fn(t) = t < t_load ? T_load :
           t < t_rest_end ? T_rest :
           min(T_rest + ramp_rate*(t - t_rest_end), 1200.0)

model = build_palioxis_trapping_model(
    palioxis_model = pal,
    mesh           = mesh,
    defects        = defects,
    temperature    = FunctionTemperature(T_fn),
    left_bc        = t -> t < t_load ? 1.0 : 0.0,
    right_bc       = t -> 0.0,
)

u0 = zeros(nvariables(model.layout) * 200)

config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    abstol = 1e-10, reltol = 1e-8,
    saveat = sort(unique(vcat(
        collect(range(0.0, t_load; step=1800.0)),
        collect(range(t_rest_end, t_rest_end + 54000.0; step=10.0)),
    ))),
)

t_end  = t_rest_end + (1200.0 - T_rest) / ramp_rate
result = wrap_result(model, solve_problem(model, u0, (0.0, t_end), config), config)
write_summary_csv(result, "tds_summary.csv")
write_field_output_hdf5(result, "tds_fields.h5")
```

---

## Visualising Results

Plotting requires a Makie backend:

```julia
using Flopsy, CairoMakie

# Desorption spectrum from summary CSV
fig = plot_tds_flux("tds_summary.csv")
save("tds_spectrum.png", fig)

# Spatial snapshot at final time (from HDF5 or in-memory result)
fig = plot_spatial_snapshot("tds_fields.h5"; time_index = :last)
save("final_profile.png", fig)

# Overlay 10 snapshots coloured by time
fig = plot_spatial_evolution(result)
save("evolution.png", fig)

# Video of the full spatial evolution
record_spatial_video(result, "spatial_evolution.mp4"; fps = 24)
```

---

## Documentation

Full documentation is built with Documenter.jl:

```
julia --project=docs docs/make.jl
```

Key pages:
- **Architecture** — state vector layout, operators, Jacobian sparsity
- **Formulations** — all formulations with worked examples
- **Hotgates Interface** — Palioxis integration, equilibrium ICs, boundary conditions
- **HDF5 Output** — file schema, Julia and Python reading examples
- **Plotting** — visualisation helpers reference
- **API Reference** — full docstring index
