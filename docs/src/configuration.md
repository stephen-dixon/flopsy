# Configuration Layer

Flopsy now treats TOML as a registered object-building language rather than a large hard-coded parser. The config system is declarative, but the top-level structure is fixed and every supported block type must be registered.

## Fixed top-level domains

Supported top-level TOML domains are:

- `mesh`
- `backend`
- `ic`
- `bc`
- `output`
- `problem`

Plugins can register new types inside those domains, but they do not invent arbitrary new top-level sections.

## Registry-backed syntax

Each syntax entry is represented by metadata in the registry:

- `domain`
- `type_name`
- parameter schema
- defaults
- required fields
- validation callback
- help text
- build callback

Core syntax and plugin syntax are registered the same way. You can inspect the live registry with:

```bash
flopsy syntax list
flopsy syntax show bc dirichlet
```

## Two-stage build process

### Stage 1: object construction

`parse_input_deck(path)` reads the TOML file into named blocks. `build_context(deck)` validates those blocks against the registry and constructs named runtime objects into a `BuildContext`.

Object groups currently include:

- meshes
- backends
- initial conditions
- boundary conditions
- outputs
- problems

### Stage 2: problem assembly

`build_simulation(ctx, problem_name)` resolves the references from `[problem.<name>]`, validates cross-object compatibility, assembles the `SystemModel`, applies ICs, and returns a ready-to-solve configuration.

For the normal user path:

```julia
using Flopsy

ctx = validate_input_deck("examples/trapping_1d_with_bc_ic.toml")
result = run_input_deck("examples/trapping_1d_with_bc_ic.toml")
```

## Example deck

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

[bc.left_H]
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
type = "trapping_1d"
mesh = "main"
backend = "main"
ics = ["mobile"]
bcs = ["left_H"]
outputs = ["fields"]
tspan = [0.0, 1.0]
saveat = [0.0, 0.25, 0.5, 0.75, 1.0]
```

## Backends define species

Flopsy does not have user-declared `[variables]` blocks yet.

Instead, each backend exposes species metadata through `SpeciesInfo`, including:

- species name
- role
- transport classification
- state index
- whether the species can be an explicit BC target

This metadata is available early enough to validate ICs and BCs before the solver is built.

## Initial conditions

Two IC styles are supported conceptually:

- species-targeting ICs such as `uniform_species`
- backend-derived ICs such as `palioxis_equilibrium`

Every IC builder must report the species it affects. Overlapping IC assignments are rejected by default.

Example species-targeting IC:

```toml
[ic.mobile]
type = "uniform_species"
species = "H_mobile"
value = 1e-3
```

Example backend-derived IC:

```toml
[ic.eq]
type = "palioxis_equilibrium"
backend = "main"
driving_quantity = "H_total"
value = 1e-3
temperature = 500.0
```

`palioxis_equilibrium` is registered from the Palioxis extension layer, not from the core built-in syntax module.

## Boundary conditions

Explicit TOML BCs must always target a named backend species:

```toml
[bc.left_H]
type = "dirichlet"
species = "H_mobile"
boundary = "left"
value = 0.0
method = "eliminated"
```

Current status:

- explicit TOML BC syntax currently supports only `dirichlet`
- if no BC blocks are supplied, the diffusion operator provides the default closed zero-flux Neumann behaviour
- BC validation rejects unknown species and stationary/internal species that do not support boundary treatment

## Outputs

The built-in output syntax is `hdf5`, with optional XDMF companion generation:

```toml
[output.fields]
type = "hdf5"
file = "fields.h5"
xdmf = true
```

You can also generate the XDMF companion later with:

```bash
flopsy xdmf fields.h5
```

## CLI workflow

Typical CLI usage:

```bash
flopsy validate input.toml
flopsy run input.toml
flopsy syntax list
flopsy syntax show mesh uniform_1d
flopsy plugin list
```

`validate` parses, validates, and assembles the named object graph without running the solver.

## Plugins

Plugins are loaded from a managed runtime Julia environment in the user config area.

This means:

- the compiled CLI does not need to be rebuilt when a plugin is installed
- plugins can register additional syntax into the same registry as built-ins
- failures are surfaced through `flopsy plugin list`

Palioxis support currently remains in this repository via `ext/PalioxisExt.jl`, but architecturally it behaves like a plugin.

## Biased mesh status

`biased_1d` is registered so the syntax is discoverable, but nonuniform mesh numerics are not implemented yet. Attempting to build it errors explicitly rather than silently using uniform formulas.
