# Input Deck Format

## Canonical Model

The supported TOML model is registry-driven and block-oriented.

Each input deck is composed of named blocks grouped under fixed top-level domains:

- `mesh`
- `backend`
- `ic`
- `bc`
- `output`
- `temperature` *(optional — required for TDS workflows)*
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

## Temperature Providers

Temperature models are declared as named `[temperature.<name>]` blocks and referenced from
`[problem.<name>]` via the `temperature` field.

```toml
[temperature.ramp]
type = "linear_ramp"
T0 = 300.0
rate = 0.5          # K/s

[problem.tds_run]
type = "tds"
mesh = "main"
backend = "main"
temperature = "ramp"
tspan = [0.0, 700.0]
dt = 1.0
```

Built-in temperature types:

| type | description |
|------|-------------|
| `constant` | Uniform temperature (`value` in K) |
| `linear_ramp` | T(t) = T0 + rate × t (`T0` K, `rate` K/s) |
| `piecewise` | Hold/ramp stages (see `stages` array) |

## Output Types

| type | description |
|------|-------------|
| `hdf5` | Full pointwise field output (`file`, optional `xdmf`) |
| `summary_csv` | Time-series scalars: integrals, surface fluxes (`file`, optional `fields` list) |

Example combining both:

```toml
[output.fields]
type = "hdf5"
file = "result.h5"
xdmf = true          # or "result.xdmf" for explicit path
export_Deff = true  # optional Palioxis effective-diffusion auxiliary field
export_equilibrium_trapped = true
export_retention_total = true
export_retention_by_occupation = true

[output.summary]
type = "summary_csv"
file = "summary.csv"
```

The equilibrium auxiliary export flags are ignored for ordinary models.  For
`backend.type = "palioxis_effective_diffusion"` they add an `/equilibrium_aux`
HDF5 group and pointwise field slices for ParaView/XDMF.

## TDS Problem Class

`type = "tds"` is a dedicated problem class for thermal desorption spectroscopy workflows.
It requires a `temperature` reference and defaults to `formulation = "split"`.

```toml
[problem.desorption]
type = "tds"
mesh = "main"
backend = "main"
temperature = "ramp"
ics = ["loaded"]
bcs = ["vacuum_left", "vacuum_right"]
outputs = ["fields", "summary"]
tspan = [0.0, 700.0]
formulation = "split"
split_method = "strang"   # or "lie"
dt = 1.0
saveat = ...
```

## Solver Options

| field | description |
|-------|-------------|
| `formulation` | `unsplit`, `imex`, `imex_reaction`, `split`, `residual` |
| `split_method` | `strang` (default) or `lie` — only used when `formulation = "split"` |
| `algorithm` | `Rodas5`, `Rodas5P`, `KenCarp4`, `CVODE_BDF` |
| `dt` | Required for `split` formulation (macro step size) |
| `abstol` | Absolute tolerance (default 1e-8) |
| `reltol` | Relative tolerance (default 1e-6) |
| `show_progress` | Show a solve progress bar (default `true`; set `false` for batch logs/tests) |

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
