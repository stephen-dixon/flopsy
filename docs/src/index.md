# Flopsy.jl Documentation

## What is Flopsy?

Flopsy is a Julia package for stiff one-dimensional reaction-diffusion problems, with a
programmable core solver framework and a thin config-driven problem layer for TOML input
decks.

The core layer is designed for direct Julia use: compose models, operators, boundary
conditions, formulations, and solver settings explicitly. The config layer parses and
validates a restricted TOML schema, selects a built-in problem template, and lowers it
into the same core APIs.

## Two-Layer Architecture

- **Core framework**: generic, type-stable Julia APIs for `SystemModel`, operators,
  formulations, and solver integration.
- **Config layer**: typed TOML parsing, validation, problem-template selection, and a
  thin execution wrapper around the core.

This keeps performance-critical code free of config parsing while still supporting
input-deck workflows for supported simulations.

## Typical Pipelines

### Programmatic core

```julia
using Flopsy
using OrdinaryDiffEq
using ADTypes: AutoFiniteDiff

model = build_trapping_model(
    mesh = Mesh1D(1.0, 101),
    k_trap = 5.0,
    k_detrap = 0.5,
    diffusion_coefficient = 0.01,
)

u0 = zeros(nvariables(model.layout) * model.context.nx)
cfg = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm = Rodas5(autodiff = AutoFiniteDiff()),
    saveat = [0.0, 0.5, 1.0],
)

sol = solve_problem(model, u0, (0.0, 1.0), cfg)
result = wrap_result(model, sol, cfg)
```

### Config-driven problem

```julia
using Flopsy

cfg = load_config("examples/trapping_1d.toml")
problem = build_problem(cfg)
result = solve(problem)
```

## Target Problems

Flopsy targets stiff, spatially-distributed hydrogen transport problems where:

- **Stiffness** arises from disparate timescales between mobile hydrogen diffusion and multi-occupancy trapping dynamics.
- **Heterogeneity** stems from spatial variation in material properties, defect densities, and temperatures.
- **Multi-occupancy trapping** requires faithful representation of defect interactions, accessible through the Palioxis library.
- **1D geometry** is sufficient for many practical materials-science applications (thin films, implantation, permeation).

## Quick Links

- **[Configuration](configuration.md)** — Typed TOML parsing, validation, and factory-built problems.
- **[Architecture](architecture.md)** — Understand the code structure: variable layouts, operators, formulations, and the solve pipeline.
- **[Hotgates Interface](hotgates.md)** — Learn how Flopsy bridges to Palioxis for accurate reaction-rate calculations.
- **API Reference** — Full API documentation (generated from docstrings).
