# Hotgates Interface

## Motivation

Palioxis provides high-fidelity multi-occupancy hydrogen trapping kinetics. Flopsy integrates that functionality without baking Palioxis-specific TOML parsing into the core registry.

Today the Palioxis implementation stays in this repository via `ext/PalioxisExt.jl`, but it is structured like a plugin:

- it registers backend syntax into the shared syntax registry
- it registers backend-derived IC syntax such as `palioxis_equilibrium`
- it exposes backend-owned species metadata to the core validator

That separation keeps Flopsy core generic while preserving Palioxis support.

## HotgatesTrappingAdaptor

The `HotgatesTrappingAdaptor` holds the configuration and metadata required to interface with Palioxis or similar trapping backends:

- **mobile_indices** — indices of mobile species in the Flopsy state vector
- **trap_indices** — indices of trapped species in the Flopsy state vector
- **mobile_names** — species names for mobile variables (metadata and postprocessing)
- **trap_names** — species names for trapped variables
- **defect_names** — species names for defect sites (not in the state vector)
- **defects** — a `(n_defects, nx)` matrix of static, spatially-varying defect densities. Each column `ix` contains the densities of all defect types at spatial node `ix`.
- **trap_groups** — `Vector{Vector{Int}}` grouping trap variable indices by defect type.
  Each inner vector lists the indices (within `trap_indices`) that belong to one defect
  type.  Used by `jacobian_node_sparsity` to construct a selective sparse Jacobian
  prototype: trap–trap entries are tridiagonal *within* each group, and dense between
  groups.  Defaults to singleton groups (each trap in its own group, fully block-diagonal).

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

`PalioxisExt` is a Julia package extension loaded when both `Flopsy` and `Palioxis` are available. It dispatches `hotgates_rates!` to the actual Palioxis library and registers Palioxis-specific TOML syntax.

### Activating PalioxisExt

To use Palioxis with Flopsy:

```julia
using Flopsy
using Palioxis
```

The extension is triggered automatically by importing `Palioxis` and makes Palioxis syntax available to the registry.

## Registered Palioxis syntax

Palioxis syntax is not hard-coded in the Flopsy core built-ins. It is registered from the extension layer.

Current extension-provided syntax includes:

- `backend.palioxis_trapping`
- `ic.palioxis_equilibrium`

If Palioxis is unavailable, those syntax entries are unavailable as well. The registry will then fail cleanly with an unknown-syntax error rather than failing later inside the solver.

You can inspect the live registry with:

```bash
flopsy syntax list
flopsy syntax show backend palioxis_trapping
flopsy syntax show ic palioxis_equilibrium
```

### Building a Palioxis Trapping Model

The extension exposes `build_palioxis_trapping_model`, which constructs a complete Flopsy
`SystemModel` backed by a `Palioxis.MultipleDefectModel`.  Diffusion coefficients are
queried automatically from Palioxis at every time step via
`Palioxis.diffusion_constants(model, T)` — no separate `D` array is needed.

```julia
using Flopsy
using Palioxis

# Create a Palioxis model with defect parameters
palioxis_model = Palioxis.MultipleDefectModel(...)

# Define spatial mesh: domain length and number of nodes
mesh = Mesh1D(1e-3, 100)   # 1 mm domain, 100 nodes

# Spatially-varying defect profile (n_traps × nx)
defects = fill(1e-3, palioxis_model.n_traps, 100)  # uniform ρ = 1e-3

# Temperature provider
T_provider = LinearRampTemperature(300.0, 10.0 / 60.0)  # 300 K, 10 K/min ramp

# Build the model — D(T) is sourced from Palioxis automatically
model = build_palioxis_trapping_model(
    palioxis_model = palioxis_model,
    mesh           = mesh,
    defects        = defects,
    temperature    = T_provider,
    left_bc        = t -> 0.0,   # vacuum at left surface
    right_bc       = t -> 0.0,   # vacuum at right surface
)
```

The model is immediately ready to use with Flopsy's solver infrastructure.

### Further Palioxis Capabilities

The Palioxis library offers capabilities beyond basic rate computations, all accessible
through the `palioxis_model` handle:

- **Temperature-dependent D** — `PalioxisDiffusionCoefficients` (created internally by
  `build_palioxis_trapping_model`) calls `Palioxis.diffusion_constants(model, T)` at every
  time step.  D is always consistent with the current temperature and model parameters.
- **Retention integrals** — `Palioxis.retention` and `Palioxis.noneq_retention` compute
  cumulative hydrogen uptake from a constant flux boundary condition.
- **Equilibrium initial conditions** — see the next section.
- **Temperature sensitivity** — Palioxis uses Arrhenius factors internally; any temperature
  profile passed through the `AbstractTemperatureProvider` interface is evaluated per time
  step and per node.

## Equilibrium Initial Conditions

For TDS desorption experiments, it is common to start from an equilibrium trapped distribution. Flopsy supports that in two ways:

- direct Julia helpers such as `build_equilibrium_ic` and `build_ic_from_total_hydrogen`
- a registered backend-derived IC syntax, `palioxis_equilibrium`

The config-driven form is useful because it can populate multiple backend-defined species without requiring user-declared variable blocks.

Example:

```toml
[backend.main]
type = "palioxis_trapping"

[ic.eq]
type = "palioxis_equilibrium"
backend = "main"
driving_quantity = "H_total"
value = 1e-3
temperature = 500.0
```

This IC validates against backend species metadata and reports the species it affects, so overlap detection works consistently with species-targeted ICs.

The underlying Julia helper path remains available:

```julia
using Flopsy
using Palioxis

# Single mobile species, length-nx vector
mobile_profile = zeros(nx)  # or a prescribed profile

# Compute trapped equilibrium at temperature T
u0 = build_equilibrium_ic(palioxis_model, model, mobile_profile, 300.0)
```

For multiple mobile species, pass a `(n_gas, nx)` matrix instead.

The function calls `Palioxis.set_initial_conditions(model, mobile_node, T)` per spatial
node to compute the equilibrium trapped concentrations given the mobile concentration
at that node.

### Equilibrium from Total Hydrogen

If you have a total-hydrogen profile (mobile + trapped summed) rather than a mobile
profile, use `build_ic_from_total_hydrogen`:

```julia
# Gaussian total-H profile centred at mid-domain
x     = model.context.mesh.x
total = 0.01 .* exp.(-(x .- x[end]/2).^2 ./ (2 * (x[end]/10)^2))

u0 = build_ic_from_total_hydrogen(palioxis_model, model, total, 400.0)
```

Per node this distributes `total[ix]` as an even initial guess, applies
`Palioxis.ensure_bounds!`, then converges to the nearest equilibrium via
`Palioxis.calculate_steady_state!`.

!!! note
    `calculate_steady_state!` finds the nearest equilibrium but does not strictly
    conserve total hydrogen.  For typical TDS conditions (strong trapping at low T)
    the result is physically accurate; verify with a mass-balance check if precision
    is critical.

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
T_provider = LinearRampTemperature(300.0, 1.0)   # T0, ramp_rate
```

Typical TDS ramp rates are 0.5–5 K/s (0.0083–0.083 K/s).

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

This is suitable for sealed boundaries or materials where external diffusion is negligible. In the TOML layer, this is the behaviour you get by omitting `[bc.*]` blocks entirely.

### Dirichlet Boundary Conditions

For surface absorption/desorption, use `WeakDirichletBoundaryOperator` or `DirichletBoundaryOperator`.

In the TOML layer, current explicit BC syntax is species-targeted `dirichlet` only. Example:

```toml
[bc.left_H]
type = "dirichlet"
species = "H_mobile"
boundary = "left"
value = 0.0
method = "weak"
```

BC validation checks that:

- the named species exists on the selected backend
- the species is a valid boundary target
- the species supports diffusive transport

```julia
boundary = WeakDirichletBoundaryOperator(
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

For stronger enforcement (e.g. when the ghost-node approach is insufficient), use
`DirichletBoundaryOperator` with `PenaltyMethod`, `MassMatrixMethod`, `CallbackMethod`,
or `EliminatedMethod`.  See **[Architecture → Boundary Condition Operators](architecture.md)**
for a comparison.

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

## Mass Conservation Diagnostics

`check_mass_conservation` verifies that the total hydrogen `H(t)` plus the cumulative
outward diffusive flux equals the initial hydrogen `H(0)`:

```julia
mc = check_mass_conservation(result; rtol=1e-3)
mc.conserved           # Bool — true if max relative error < rtol
mc.max_relative_error  # Float64 — worst-case relative balance error
mc.total_hydrogen      # Vector{Float64} — H(t) at each saved time
mc.cumulative_flux     # Vector{Float64} — ∫ flux dτ from t[1] to t[k]
mc.balance             # Vector{Float64} — H + cumulative_flux (should ≈ H[1])
```

The flux uses the mass-conserving formula `D*(U[boundary] − g)/dx` from each active
`WeakDirichletBoundaryOperator`.  Closed systems (Neumann BCs) have no boundary operator
and the cumulative flux is identically zero; `mc.conserved` is `true` if H is constant
to within `rtol`.

!!! note "Temporal resolution"
    `check_mass_conservation` integrates boundary fluxes by the trapezoidal rule over the
    saved time points.  Accurate mass balance requires that `saveat` is fine enough to
    resolve the diffusion dynamics (`dt_save ≪ L²/(π²D)`).  For fast-diffusing species
    or coarse `saveat`, the trapezoidal error can exceed `rtol`.

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
u0_implant = ...
tspan_implant = ...
config_implant = ...
sol_implant = solve_problem(model_implant, u0_implant, tspan_implant, config_implant)
result_implant = wrap_result(model_implant, sol_implant, config_implant)
write_field_output_hdf5(result_implant, "implant_out.h5")

# Stage 2: TDS from implanted state
model_tds = ...  # may have different variable set or mesh
config_tds = ...
tspan_tds = ...
u0 = load_ic_from_hdf5("implant_out.h5", model_tds)
sol_tds = solve_problem(model_tds, u0, tspan_tds, config_tds)
result_tds = wrap_result(model_tds, sol_tds, config_tds)
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
