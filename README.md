# flopsy

![Logo](./docs/flopsy.png)

Flopsy is a Julia package for stiff 1D reaction-diffusion simulations with two user paths:

- a programmable core API for building models directly in Julia
- a registry-driven TOML layer and `flopsy` CLI for supported input-deck workflows

The config path is a thin wrapper over the same solver stack. It does not define a second execution architecture.

## Design

Flopsy now builds input decks through a registered syntax system inspired by MOOSE.
Built-ins and plugins register syntax entries into the same in-memory registry, keyed by fixed domain and type pairs such as:

- `(:mesh, :uniform_1d)`
- `(:backend, :trapping_1d)`
- `(:ic, :uniform_species)`
- `(:ic, :palioxis_equilibrium)`
- `(:bc, :dirichlet)`
- `(:output, :hdf5)`
- `(:problem, :simulation)`

The top-level TOML structure is fixed:

- `[mesh.<name>]`
- `[backend.<name>]`
- `[ic.<name>]`
- `[bc.<name>]`
- `[output.<name>]`
- `[problem.<name>]`

Backends own the species/state layout. Users do not declare `[variables]` blocks yet. Initial conditions and explicit boundary conditions must target backend-declared species by name.

## Quickstart

### Programmatic core

```julia
using Flopsy
using OrdinaryDiffEq

mesh = Mesh1D(1e-3, 60)
model = build_trapping_model(
    mesh = mesh,
    k_trap = 5e-3,
    k_detrap = 0.05,
    diffusion_coefficient = 1e-7,
)

u0 = zeros(nvariables(model.layout) * length(mesh.x))
tspan = (0.0, 200.0)

config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm = Rodas5P(),
    abstol = 1e-9,
    reltol = 1e-7,
    saveat = collect(0.0:2.0:200.0),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
write_field_output_hdf5(result, "fields.h5")
write_xdmf_for_flopsy_h5("fields.h5")
```

### Registry-driven TOML

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

Run it with:

```bash
julia --project=. scripts/run_from_toml.jl examples/trapping_1d_with_bc_ic.toml
```

or via the CLI entrypoint:

```bash
flopsy run examples/trapping_1d_with_bc_ic.toml
```

## Boundary conditions and initial conditions

Current TOML BC support is intentionally narrow:

- explicit BC syntax currently supports `dirichlet`
- if no `[bc.*]` blocks are provided, diffusion uses the default closed zero-flux Neumann behaviour

Examples:

- [examples/trapping_1d_with_bc_ic.toml](/Users/sdixon/src/palioxis-tds/flopsy/examples/trapping_1d_with_bc_ic.toml)
- [examples/closed_system_no_bc.toml](/Users/sdixon/src/palioxis-tds/flopsy/examples/closed_system_no_bc.toml)

ICs come in two categories:

- species-targeted ICs such as `uniform_species`
- backend-derived ICs such as `palioxis_equilibrium`

Overlapping IC assignments are rejected by default.

## CLI

The CLI is designed to work from source and as a PackageCompiler target.

```text
flopsy run <input.toml>
flopsy validate <input.toml>
flopsy syntax list
flopsy syntax show <domain> <type>
flopsy plugin list
flopsy plugin register <name> --registry <url>
flopsy plugin register <name> --url <pkg-url>
flopsy plugin register <name> --path <local-path>
flopsy plugin remove <name>
flopsy xdmf <fields.h5> [--output out.xdmf]
```

`validate` parses and builds the named objects without running the solver. `syntax list` and `syntax show` are generated from the live registry, so built-ins and loaded plugins appear through the same interface.

## Plugins

Flopsy supports runtime-installed plugins through a managed Julia environment under the user config area. The compiled CLI does not mutate or rebuild itself when plugins are installed.

At startup Flopsy:

- builds the built-in syntax registry
- activates and inspects the managed plugin environment
- loads plugin packages
- lets each plugin call `register_flopsy_plugin!(registry)`

Palioxis still lives in this repository today through `ext/PalioxisExt.jl`, but it is wired through the same registration model as a plugin.

## XDMF and ParaView

HDF5 field outputs can request an XDMF companion directly:

```toml
[output.fields]
type = "hdf5"
file = "fields.h5"
xdmf = true
```

You can also generate one afterwards:

```bash
flopsy xdmf fields.h5
```

This writes a `.xdmf` file that points ParaView at the HDF5 field datasets.

## Current limitations

- explicit TOML BC syntax is `dirichlet` only
- omitted BC blocks are the supported path for closed zero-flux boundaries
- `biased_1d` mesh syntax is registered, but nonuniform numerics are not implemented yet and it errors explicitly
- the old typed `load_config` and `build_problem(::ProblemConfig)` path still exists for compatibility, but the registry-driven path is now the primary public interface

## Documentation

Build the docs with:

```bash
julia --project=docs docs/make.jl
```

Key pages:

- [Configuration](docs/src/configuration.md)
- [Architecture](docs/src/architecture.md)
- [Hotgates Interface](docs/src/hotgates.md)
- [HDF5 Output](docs/src/hdf5.md)
- [API Reference](docs/src/api.md)
