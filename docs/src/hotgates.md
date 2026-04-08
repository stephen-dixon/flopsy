# Hotgates Interface

## Motivation

Palioxis is a high-fidelity Fortran library for modeling multi-occupancy hydrogen trapping in materials. It accurately captures the complex kinetics of hydrogen binding to defect sites, detrapping barriers, and site saturation effects. The Hotgates library provides language bindings that make Palioxis accessible from higher-level environments.

Flopsy integrates Palioxis trapping kinetics through a dedicated adapter layer: the `HotgatesReactionOperator` and `HotgatesTrappingAdaptor`. This integration allows spatially-resolved 1D reaction-diffusion simulations to benefit from Palioxis's sophisticated defect models while remaining composable with other Flopsy operators.

## HotgatesTrappingAdaptor

The `HotgatesTrappingAdaptor` holds the configuration and metadata required to interface with Palioxis or similar trapping backends:

- **mobile_indices** — indices of mobile species in the Flopsy state vector
- **trap_indices** — indices of trapped species in the Flopsy state vector
- **mobile_names** — species names for mobile variables (metadata and postprocessing)
- **trap_names** — species names for trapped variables
- **defect_names** — species names for defect sites (not in the state vector)
- **defects** — a `(n_defects, nx)` matrix of static, spatially-varying defect densities. Each column `ix` contains the densities of all defect types at spatial node `ix`.

The key design choice is that **defect densities are spatially varying but static in time**. This allows the trapping model to be computed independently at each node without global communication, while still respecting material heterogeneity.

## HotgatesReactionOperator

The `HotgatesReactionOperator` is a reaction operator that computes node-local trapping and detrapping rates using the Hotgates interface.

For each spatial node `ix`:

1. Extract the mobile hydrogen concentrations at that node.
2. Extract the trapped hydrogen concentrations (one per defect type) at that node.
3. Extract the trap densities from the adaptor's `defects` matrix, column `ix`.
4. Query the current temperature via the temperature provider.
5. Call `hotgates_rates!(dmobile, dtrapped, model, mobile, defects, trapped, T)` to compute the rate of change.
6. Accumulate the computed rates into the full `du` vector.

Because this computation is node-local, the operator scales well to fine spatial discretizations and can be parallelized in the future.

## hotgates_rates! Dispatch

The function `hotgates_rates!` is the key extension point for backend integration. Its signature is:

```julia
hotgates_rates!(dmobile, dtrapped, model, mobile, defects, trapped, T)
```

where:
- `dmobile` — output vector of rates of change for each mobile species (written, not accumulated).
- `dtrapped` — output vector of rates for each trapped species (written, not accumulated).
- `model` — the Palioxis model instance, or a mock for testing.
- `mobile` — current mobile hydrogen concentrations at this node.
- `defects` — trap densities at this node.
- `trapped` — trapped hydrogen concentrations at this node.
- `T` — temperature at this node in Kelvin.

### FakeHotgatesModel

For testing and validation without a full Palioxis installation, Flopsy provides `FakeHotgatesModel`. It implements simple mass-action kinetics with configurable rate constants:

```julia
model = FakeHotgatesModel(k_trap=1e-3, k_detrap=0.1)
```

The rates are computed as:

```
trap_rate   = k_trap * mobile * (1 - trapped)
detrap_rate = k_detrap * trapped
```

This mock allows unit tests and diagnostic runs to proceed without external dependencies.

## PalioxisExt

The `PalioxisExt` is a Julia 1.9+ package extension (loaded automatically when both `Palioxis` and `Flopsy` are imported). It dispatches `hotgates_rates!` to the actual Palioxis library.

### Activating PalioxisExt

To use Palioxis with Flopsy:

```julia
using Flopsy
using Palioxis
```

The extension is triggered automatically by importing `Palioxis` and connects Flopsy to the Palioxis backend.

### Building a Palioxis Trapping Model

The extension exposes `build_palioxis_trapping_model`, which constructs a complete Flopsy `SystemModel` backed by a `Palioxis.MultipleDefectModel`:

```julia
using Flopsy
using Palioxis

# Create a Palioxis model with defect parameters
palioxis_model = Palioxis.MultipleDefectModel(...)

# Define spatial mesh
mesh = Mesh1D(0.0, 1e-3, 100)

# Spatially-varying defect profile (n_traps × nx)
defects = ...  # build using Palioxis.get_defect_concentrations per node

# Temperature provider
T_provider = LinearRampTemperature(T0=300.0, ramp_rate=1.0)

# Diffusion coefficients (one per mobile species, zeros for traps)
D = [1e-8, 1e-8]  # m²/s

# Build the model
model = build_palioxis_trapping_model(
    palioxis_model=palioxis_model,
    mesh=mesh,
    defects=defects,
    temperature=T_provider,
    diffusion_coefficients=D,
)
```

The model is immediately ready to use with Flopsy's solver infrastructure.

### Further Palioxis Capabilities

The Palioxis library offers capabilities beyond basic rate computations, all accessible through the `palioxis_model` handle:

- **Jacobians** — `Palioxis.time_derivatives_jacobian` computes sensitivities of rates with respect to mobile and trapped concentrations. This enables efficient implicit solving.
- **Retention integrals** — `Palioxis.retention` and `Palioxis.noneq_retention` compute the cumulative hydrogen uptake (moles per unit area) from a constant flux boundary condition. Useful for comparison with permeation experiments.
- **Equilibrium initial conditions** — see the next section.
- **Temperature sensitivity** — Palioxis models temperature through standard Arrhenius factors. Any temperature profile can be passed at each time step, allowing thermal transients to be captured.

## Equilibrium Initial Conditions

For TDS desorption experiments, it is common to start with hydrogen trapped in defects at some initial temperature, having negligible mobile concentration. The `build_equilibrium_ic` function sets trapped species to Palioxis equilibrium given a mobile profile:

```julia
using Flopsy
using Palioxis

# Single mobile species, length-nx vector
mobile_profile = zeros(nx)  # or a prescribed profile

# Compute trapped equilibrium at temperature T
u0 = build_equilibrium_ic(
    palioxis_model,
    adaptor,
    mobile_profile,
    T=300.0,
)
```

For multiple mobile species, pass a `(n_gas, nx)` matrix instead.

The function calls `Palioxis.set_initial_conditions(model, mobile_node, T)` per spatial node to compute the equilibrium trapped concentrations given the mobile concentration at that node.

### Achieving Full Equilibrium

If you need the **mobile and trapped species both at equilibrium** simultaneously (e.g., to conserve total hydrogen during a cold start), use `Palioxis.calculate_steady_state!` on each node manually and construct the state vector yourself.

## Temperature Providers

Temperature is supplied via an `AbstractTemperatureProvider`. Three implementations are available.

### ConstantTemperature

Returns the same temperature everywhere:

```julia
T_provider = ConstantTemperature(300.0)  # K
```

### LinearRampTemperature

Implements a linear temperature ramp, standard for TDS experiments:

```julia
T(t) = T0 + ramp_rate * t
```

```julia
# Ramp from 300 K at 1 K/s
T_provider = LinearRampTemperature(T0=300.0, ramp_rate=1.0)
```

Typical TDS ramp rates are 0.5–5 K/s.

### FunctionTemperature

For arbitrary callable `f(t)` returning Kelvin:

```julia
# Ramp to 1000 K over 700 s, then hold
T_func(t) = min(300.0 + t, 1000.0)
T_provider = FunctionTemperature(T_func)
```

The same temperature value is used at every spatial node (uniform temperature field across the domain).

## Boundary Conditions

### Neumann (Zero-Flux) Boundary Conditions

By default, `LinearDiffusionOperator` implements zero-flux (Neumann) boundary conditions at both surfaces:

```
∂u/∂n = 0
```

This is suitable for sealed boundaries or materials where external diffusion is negligible.

### Dirichlet Boundary Conditions

For surface absorption/desorption, use `DirichletBoundaryOperator` to impose fixed surface concentrations:

```julia
boundary = DirichletBoundaryOperator(
    selector,
    diffusion_coefficients;
    left  = t -> 0.0,        # vacuum at left surface
    right = t -> 0.0,        # vacuum at right surface
)
```

The operator adds ghost-node corrections on top of the Neumann stencil already in `LinearDiffusionOperator`. The combined stencil at a boundary becomes:

```
dU[1] += D * (U[2] - 2*U[1] + g) / dx²
```

where `g = left(t)` or `g = right(t)` is the Dirichlet value.

### Time-Varying Boundary Conditions

The boundary callables receive time `t` and can implement time-varying conditions. This is useful for combined implantation and desorption in a single run:

```julia
# Implantation phase until t_cutover, then vacuum
boundary = DirichletBoundaryOperator(
    selector,
    diffusion_coefficients;
    left = t -> t < t_cutover ? source_conc : 0.0,
    right = t -> 0.0,
)
```

!!! note
    Set the initial condition at boundary nodes to match `left(0)` and `right(0)` to avoid discontinuities at t=0.

### Alternative: Chaining Simulations

For complex multi-phase experiments, the alternative to time-varying boundary conditions is to run separate simulations and chain them via HDF5 (see the next sections). For example: implantation with a fixed source, followed by desorption from the implanted state.

## Surface Flux Summaries

The function `surface_diffusive_fluxes` queries the `LinearDiffusionOperator` and computes the outward diffusive flux at both boundary nodes for each diffusing variable over all saved times.

```julia
fluxes = surface_diffusive_fluxes(result)
# fluxes::Dict{Symbol, NamedTuple}
# fluxes[:mobile_H] = (left=Vector, right=Vector)
```

### Sign Convention

Fluxes are **positive in the direction of the outward-facing normal**:

- **Left boundary** (x=0):  `left_flux = D * (C[2] - C[1]) / dx`
  - Normal points in −x direction
  - Equals `D * ∂C/∂x|_{x=0}`
  - Positive during TDS desorption (H leaving the sample at left surface)

- **Right boundary** (x=L):  `right_flux = D * (C[nx-1] - C[nx]) / dx`
  - Normal points in +x direction
  - Equals `−D * ∂C/∂x|_{x=L}`
  - Positive during TDS desorption (H leaving the sample at right surface)

### Automatic CSV Output

The fluxes are automatically computed and included in the summary CSV generated by `build_summary_dataframe` and `write_summary_csv`, with columns:
- `left_flux_<varname>`
- `right_flux_<varname>`

This is convenient for postprocessing desorption experiments.

## Chaining Simulations via HDF5

For multi-stage experiments (e.g., implantation → resting period → TDS), Flopsy supports loading initial conditions from previous simulation output via `load_ic_from_hdf5`.

### Overview

1. **Stage 1** — Run implantation simulation, save the final state to HDF5
2. **Stage 2** — Load the HDF5 output as the initial condition for a subsequent simulation

### Workflow

```julia
using Flopsy

# Stage 1: Implantation
model_implant = ...
config_implant = ...
result_implant = run_simulation(model_implant, config_implant)
write_field_output_hdf5(result_implant, "implant_out.h5")

# Stage 2: TDS from implanted state
model_tds = ...  # may have different variable set or mesh
u0 = load_ic_from_hdf5("implant_out.h5", model_tds)
# u0 is ready to pass to solve_problem
result_tds = solve_problem(
    build_problem(model_tds, config_tds, u0),
    config_tds,
)
```

### Behavior

`load_ic_from_hdf5(path, model; time_index=:last)` reads the HDF5 file written by `write_field_output_hdf5` and constructs a Flopsy state vector:

- Matches variables by name between the file and `model.layout`
- Variables present in the model but absent from the file are initialized to zero with a warning
- This allows chaining simulations where the variable set differs (e.g., adding a new trap type after implantation)
- The mesh node count must match; this is checked and raises `DimensionMismatch` if they differ

`time_index` selects which saved time to load:
- `:last` (default) — load the final state
- An integer — load state at that saved-time index

### HDF5 File Structure

Files written by `write_field_output_hdf5` have the structure:

```
/time              — saved times (Float64 vector)
/mesh/x            — spatial node coordinates
/mesh/dx           — node spacing
/fields/<varname>  — matrix (nt, nx) for each variable
/metadata/...      — selected metadata as string datasets
```
