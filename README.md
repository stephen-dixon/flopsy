# flopsy

![Logo](./docs/flopsy.png)

Flopsy is a Julia package for stiff 1D reaction-diffusion simulations.

It supports two primary use cases:

- running registry-driven TOML input decks through the CLI
- using the package directly as a Julia library/framework

The old typed TOML/config path still exists only as a deprecated compatibility layer. The supported TOML workflow is the registry-driven input-deck system built from named blocks and registered syntax.

## Installation and Development

Clone the repository and instantiate the environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

For source-based CLI use, run:

```bash
julia --project=. scripts/flopsy.jl --help
```

## Quick Start: CLI

The canonical input-deck format is the registry-driven named-block TOML model.

Minimal example: [`examples/minimal_input_deck.toml`](examples/minimal_input_deck.toml)

Validate it:

```bash
julia --project=. scripts/flopsy.jl validate examples/minimal_input_deck.toml
```

Run it:

```bash
julia --project=. scripts/flopsy.jl run examples/minimal_input_deck.toml
```

List available built-in and plugin-provided syntax:

```bash
julia --project=. scripts/flopsy.jl syntax list
```

Show one syntax entry and its schema:

```bash
julia --project=. scripts/flopsy.jl syntax show bc dirichlet
```

For a more complete deck with BC/IC/output wiring, see [`examples/trapping_1d_with_bc_ic.toml`](examples/trapping_1d_with_bc_ic.toml).

## Quick Start: Julia Library

### Run an input deck from Julia

```julia
using Flopsy

ctx = validate_input_deck("examples/minimal_input_deck.toml")
result = run_input_deck("examples/minimal_input_deck.toml")
```

Equivalent runnable example: [`examples/run_input_deck_from_julia.jl`](examples/run_input_deck_from_julia.jl)

### Build and solve directly from Julia

```julia
using Flopsy
using OrdinaryDiffEq

mesh = Mesh1D(1.0, 51)
model = build_trapping_model(
    mesh = mesh,
    k_trap = 2.0,
    k_detrap = 0.25,
    diffusion_coefficient = 0.05,
)

u0 = zeros(nvariables(model.layout) * length(mesh.x))
tspan = (0.0, 0.25)

config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm = Rodas5P(),
    abstol = 1e-8,
    reltol = 1e-6,
    saveat = [0.0, 0.125, 0.25],
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

Equivalent runnable example: [`examples/programmatic_trapping.jl`](examples/programmatic_trapping.jl)

Use the CLI when you want a stable declarative workflow, validation, and plugin-driven syntax discovery. Use the Julia API when you want programmatic model assembly, custom orchestration, or lower-level extension hooks.

## Input Deck Overview

Supported top-level TOML domains are fixed:

- `mesh`
- `backend`
- `ic`
- `bc`
- `output`
- `problem`

Each named block declares a registered `type`. Built-ins and plugins both register syntax through the same `SyntaxRegistry`.

Example:

```toml
[mesh.main]
type = "uniform_1d"
xmin = 0.0
xmax = 1.0
nx = 201

[backend.main]
type = "trapping_1d"
k_trap = 5.0
k_detrap = 0.5
diffusion_coefficient = 0.01

[ic.mobile]
type = "uniform_species"
species = "H_mobile"
value = 1e-3

[output.fields]
type = "hdf5"
file = "trapping_fields.h5"
xdmf = true

[problem.run]
type = "simulation"
mesh = "main"
backend = "main"
ics = ["mobile"]
outputs = ["fields"]
tspan = [0.0, 1.0]
saveat = [0.0, 0.5, 1.0]
```

Validation is schema-driven:

- unknown keys are rejected by default
- required fields are enforced centrally
- field types and enum-like values are checked before builders run
- object and species reference failures are surfaced with named diagnostics

## Plugins and Extensibility

Plugins extend the same syntax registry used by built-ins. A plugin package registers syntax by defining:

```julia
register_flopsy_plugin!(registry)
```

Users can manage runtime plugins through the CLI:

```bash
julia --project=. scripts/flopsy.jl plugin list
julia --project=. scripts/flopsy.jl plugin register MyPlugin --path /path/to/MyPlugin
julia --project=. scripts/flopsy.jl plugin remove MyPlugin
```

Built-in syntax and plugin syntax coexist in one live registry, so `syntax list` and `syntax show` reflect both.

## Standalone CLI Compilation

The CLI is usable from source and also prepared for PackageCompiler app builds.

App wrapper package: [`app/FlopsyCLI`](app/FlopsyCLI)

Build script:

```bash
julia --project=. scripts/build_app.jl
```

This uses `PackageCompiler.create_app` to build a standalone executable named `flopsy` under `build/flopsy-app`.

## Legacy Configuration API

The following typed-config APIs are deprecated and are no longer the recommended TOML workflow:

- `load_config`
- `parse_config`
- `validate(::ProblemConfig)`
- `build_problem(::ProblemConfig)`

They remain temporarily for compatibility, but new TOML-facing features should go through the registry/input-deck path instead.

## Documentation

Build the docs with:

```bash
julia --project=docs docs/make.jl
```

Start with:

- [`docs/src/index.md`](docs/src/index.md)
- [`docs/src/cli.md`](docs/src/cli.md)
- [`docs/src/julia_library.md`](docs/src/julia_library.md)
- [`docs/src/input_deck.md`](docs/src/input_deck.md)
- [`docs/src/plugins.md`](docs/src/plugins.md)
- [`docs/src/standalone_compilation.md`](docs/src/standalone_compilation.md)
