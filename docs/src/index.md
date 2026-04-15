# Flopsy.jl Documentation

## What is Flopsy?

Flopsy is a Julia package for stiff one-dimensional reaction-diffusion simulations with:

- a programmable, type-stable core for direct Julia use
- a registry-driven TOML layer and CLI for supported input-deck workflows

Both paths lower to the same `SystemModel` and solver APIs.

## Registered syntax model

Input decks are built from a fixed set of top-level domains:

- `mesh`
- `backend`
- `ic`
- `bc`
- `output`
- `problem`

Each named block declares a registered `type`. The registry is keyed by `(domain, type)` pairs, and each entry carries:

- a parameter schema
- defaults and required fields
- validation logic
- help text
- a build callback

Built-ins register syntax through the same mechanism as plugins. Unknown syntax fails with a registry-based error instead of falling through a hard-coded parser branch.

## Build pipeline

The config path has two stages.

### Stage 1: build named objects

Flopsy parses TOML into named blocks such as `[mesh.main]` or `[ic.eq]`, validates those blocks against the syntax registry, and stores built objects into a `BuildContext`.

### Stage 2: assemble the simulation

A `[problem.<name>]` block references the previously built mesh, backend, IC, BC, and output objects. Flopsy then assembles a `SimulationProblem` and solves it through the standard solver stack.

## Species are backend-defined

Backends own the state layout and expose species metadata through `SpeciesInfo`.

That metadata is used to validate:

- which species exist
- which species support explicit BC treatment
- which species each IC affects
- whether multiple ICs overlap

Flopsy does not expose user-declared `[variables]` TOML blocks yet. That is deliberate: backends are the source of truth for the state layout.

## Typical workflow

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
file = "fields.h5"
xdmf = true

[problem.run]
type = "trapping_1d"
mesh = "main"
backend = "main"
ics = ["mobile"]
outputs = ["fields"]
tspan = [0.0, 1.0]
saveat = [0.0, 0.5, 1.0]
```

Then:

```bash
flopsy validate input.toml
flopsy run input.toml
```

## CLI and plugins

The `flopsy` CLI supports:

- `run`
- `validate`
- `syntax list`
- `syntax show`
- `plugin list`
- `plugin register`
- `plugin remove`
- `xdmf`

Plugins are installed into a managed runtime Julia environment. This keeps the compiled CLI stable while still allowing new syntax to be added after installation.

## Current status

- explicit TOML BC syntax currently supports `dirichlet`
- if no `[bc.*]` blocks are supplied, diffusion uses implicit zero-flux Neumann behaviour
- `biased_1d` is registered, but nonuniform numerics are not implemented yet
- Palioxis equilibrium IC registration lives in `ext/PalioxisExt.jl`

## Quick links

- [Configuration](configuration.md)
- [Architecture](architecture.md)
- [Hotgates Interface](hotgates.md)
- [HDF5 Output](hdf5.md)
- [API Reference](api.md)
