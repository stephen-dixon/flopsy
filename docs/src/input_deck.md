# Input Deck Format

## Canonical Model

The supported TOML model is registry-driven and block-oriented.

Each input deck is composed of named blocks grouped under fixed top-level domains:

- `mesh`
- `backend`
- `ic`
- `bc`
- `output`
- `problem`

Each block must declare a `type` that matches a registered syntax entry in the live `SyntaxRegistry`.

## Named Blocks and References

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

[bc.left_mobile]
type = "dirichlet"
species = "H_mobile"
boundary = "left"
value = 0.0
method = "weak"

[output.fields]
type = "hdf5"
file = "trapping_fields.h5"
xdmf = true

[problem.run]
type = "simulation"
mesh = "main"
backend = "main"
ics = ["mobile"]
bcs = ["left_mobile"]
outputs = ["fields"]
tspan = [0.0, 1.0]
saveat = [0.0, 0.25, 0.5, 0.75, 1.0]
```

The `[problem.<name>]` block references previously defined objects by name.

## Validation Rules

Validation is enforced centrally from each syntax entry's `ParameterSpec` schema before builder code runs.

By default, Flopsy:

- rejects unknown keys
- enforces required keys
- checks field types where declared
- checks enum-like values where declared
- reports the failing block and field in the error message

Cross-object validation then checks:

- referenced meshes, backends, ICs, BCs, and outputs exist
- IC targets resolve to backend-defined species
- BC targets resolve to backend-defined species
- unsupported BC targets are rejected
- overlapping IC assignments are rejected

## Backend-Owned Species

Backends define the species/state layout. Users do not declare a separate `[variables]` section.

That means:

- species names are backend-defined
- IC and BC target names are validated against backend metadata
- plugin-provided backends can extend the species model without changing the top-level deck structure

## Inspecting Available Syntax

The registry is discoverable from the CLI:

```bash
julia --project=. scripts/flopsy.jl syntax list
julia --project=. scripts/flopsy.jl syntax show backend trapping_1d
```

## Validation From Julia

```julia
using Flopsy

deck = parse_input_deck("examples/minimal_input_deck.toml")
ctx = build_context(deck)
plan = build_simulation(ctx)
```

This gives direct access to the parsed deck, built object graph, and assembled simulation.
