# Configuration Layer

Flopsy now exposes a dedicated config-driven problem layer on top of the core solver
framework. The intent is:

- users can run supported simulations from TOML input decks
- developers can still build models directly in Julia
- both paths execute through the same `SystemModel` and solver machinery

## Pipeline

The config path is intentionally thin:

```julia
using Flopsy

cfg = load_config("examples/trapping_1d.toml")
problem = build_problem(cfg)
result = solve(problem)
```

`run_simulation(path)` is just a convenience wrapper over those three steps.

## Typed Config Objects

The config layer parses TOML into typed structs instead of passing raw dictionaries
through the solver stack:

- `MeshConfig`
- `BoundaryConditionConfig`
- `InputSolverConfig`
- `ProblemConfig`

This keeps parsing and validation in the outer layer while the inner solver code stays
fully typed.

## Input Schema

The built-in problem templates expect a structure like:

```toml
problem_type = "trapping_1d"

[mesh]
kind = "uniform_1d"
nx = 101
xmin = 0.0
xmax = 1.0

[solver]
formulation = "unsplit"
algorithm = "Rodas5"
abstol = 1e-8
reltol = 1e-6
saveat = [0.0, 0.25, 0.5, 0.75, 1.0]

[parameters]
t0 = 0.0
tend = 1.0
k_trap = 5.0
k_detrap = 0.5
diffusion_coefficient = 0.01
initial_mobile_pulse_amplitude = 1.0
initial_trap_occupancy = 0.0
```

Examples are provided in:

- `examples/diffusion_1d.toml`
- `examples/trapping_1d.toml`
- `examples/hotgates_trapping.toml`

## Validation

`validate(cfg)` checks the typed configuration before any solver internals are built. It
is responsible for:

- mesh sanity checks
- supported problem type checks
- solver tolerance and time-span checks
- split-formulation `dt` requirements
- boundary-condition schema checks

The goal is to fail early with user-facing errors instead of surfacing invalid input from
inside SciML or operator code.

## Problem Templates

Problem selection is handled through explicit templates, not a monolithic branch-heavy
runner:

- `Diffusion1DTemplate`
- `Trapping1DTemplate`
- `HotgatesTrappingTemplate`

The factory is responsible for choosing a template. Each template is responsible for:

- building the mesh
- assembling the `SystemModel`
- creating initial conditions
- lowering `InputSolverConfig` into core `SolverConfig`

## Legacy Input Decks

The parser still normalises the older flat TOML examples into the new typed schema where
possible. New documentation and examples use the structured `[mesh]`, `[solver]`, and
`[parameters]` layout.
