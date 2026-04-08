"""
Package extension that wires the real `Palioxis.MultipleDefectModel` into
Flopsy's Hotgates adapter.  Activated automatically when both `Flopsy` and
`Palioxis` are loaded in the same Julia session.
"""
module PalioxisExt

using Flopsy
using Palioxis

# ---------------------------------------------------------------------------
# hotgates_rates! — core dispatch for the real Palioxis backend
# ---------------------------------------------------------------------------

"""
    hotgates_rates!(dmobile, dtrapped, model::Palioxis.MultipleDefectModel,
                    mobile, defects, trapped, T)

Delegates to `Palioxis.rates!` to compute the time derivatives of mobile and
trapped species at a single spatial node.

`dmobile` and `dtrapped` are overwritten (not accumulated).
`defects` is the per-node defect concentration vector of length `model.n_traps`.
`T` is temperature in Kelvin.
"""
function Flopsy.hotgates_rates!(
    dmobile::AbstractVector,
    dtrapped::AbstractVector,
    model::Palioxis.MultipleDefectModel,
    mobile::AbstractVector,
    defects::AbstractVector,
    trapped::AbstractVector,
    T::Real,
)
    Palioxis.rates!(dmobile, dtrapped, model, mobile, defects, trapped, T)
    return nothing
end

# ---------------------------------------------------------------------------
# build_palioxis_trapping_model — convenience model builder
# ---------------------------------------------------------------------------

"""
    build_palioxis_trapping_model(; palioxis_model, mesh, defects, temperature,
                                    diffusion_coefficients, boundary=nothing)

Build a Flopsy `SystemModel` backed by a real `Palioxis.MultipleDefectModel`.

# Arguments
- `palioxis_model`         — a fully-constructed `Palioxis.MultipleDefectModel`
- `mesh`                   — `Mesh1D` defining the spatial domain
- `defects`                — `AbstractMatrix` of shape `(n_traps, nx)` holding the
                             spatially-varying defect concentration profile.
                             Column `ix` is the defect vector at mesh node `ix`.
                             Build with `Palioxis.get_defect_concentrations(model, z1, z2)`
                             for each node interval.
- `temperature`            — an `AbstractTemperatureProvider` (e.g. `ConstantTemperature`,
                             `LinearRampTemperature`, `FunctionTemperature`)
- `diffusion_coefficients` — length-`n_gas` vector of diffusion coefficients [m²/s].
                             Trapped species do not diffuse; do not include them here.
- `boundary`               — optional `DirichletBoundaryOperator` for surface BCs
                             (default `nothing` = zero-flux / Neumann at both surfaces)

# Variable layout
The resulting model has `n_gas + n_ne_species` variables per node, ordered as:
  1..n_gas          — mobile species  (tagged :mobile, :reaction, :diffusion)
  n_gas+1..end      — trapped species (tagged :trap,   :reaction)

Species names are taken directly from the Palioxis model metadata and appear in
Flopsy's `VariableLayout`, HDF5 output, and summary CSV.

# Retention / Jacobian
The underlying `Palioxis.MultipleDefectModel` exposes further capabilities
(Jacobians, retention integrals, equilibrium ICs) through the `Palioxis` API.
Pass the same `palioxis_model` handle to those functions as needed.
"""
function Flopsy.build_palioxis_trapping_model(;
    palioxis_model::Palioxis.MultipleDefectModel,
    mesh::Flopsy.Mesh1D,
    defects::AbstractMatrix,
    temperature::Flopsy.AbstractTemperatureProvider,
    diffusion_coefficients::AbstractVector,
    boundary = nothing,
)
    nx = length(mesh.x)

    size(defects) == (palioxis_model.n_traps, nx) ||
        throw(DimensionMismatch(
            "defects must have shape ($(palioxis_model.n_traps), $nx), " *
            "got $(size(defects))"
        ))

    length(diffusion_coefficients) == palioxis_model.n_gas ||
        throw(ArgumentError(
            "diffusion_coefficients must have length n_gas=$(palioxis_model.n_gas); " *
            "trapped species do not diffuse"
        ))

    n_gas = palioxis_model.n_gas
    n_ne  = palioxis_model.n_ne_species

    adaptor = Flopsy.HotgatesTrappingAdaptor(
        collect(1:n_gas),
        collect(n_gas+1 : n_gas+n_ne),
        Palioxis.gas_names(palioxis_model),
        Palioxis.trap_names(palioxis_model),
        Palioxis.defect_names(palioxis_model),
        Matrix{Float64}(defects),
    )

    # Diffusion only for mobile species; trap slots get zero
    diff_coeffs = vcat(collect(Float64, diffusion_coefficients), zeros(Float64, n_ne))

    return Flopsy.build_hotgates_trapping_model(
        mesh                   = mesh,
        model                  = palioxis_model,
        adaptor                = adaptor,
        temperature            = temperature,
        diffusion_coefficients = diff_coeffs,
        boundary               = boundary,
    )
end

# ---------------------------------------------------------------------------
# build_equilibrium_ic — initial condition from mobile profile + Palioxis equilibrium
# ---------------------------------------------------------------------------

"""
    build_equilibrium_ic(palioxis_model, adaptor, mobile_profile, T) -> Vector{Float64}

Build a flat Flopsy initial-state vector where mobile species are set from
`mobile_profile` and trapped species are placed at their Palioxis equilibrium
values for temperature `T` (in Kelvin).

`mobile_profile` may be:
- A `Vector` of length `nx` (single mobile species: `n_gas == 1`)
- A `Matrix` of shape `(n_gas, nx)` (multiple mobile species)

Uses `Palioxis.set_initial_conditions(model, mobile_node, T)` per spatial node,
which returns the equilibrium trapped concentrations given the mobile concentration.

# Notes
- For TDS initial conditions the typical workflow is: compute or prescribe the
  total implanted hydrogen profile, assume most is trapped at the starting
  temperature, pass a small mobile concentration (often ≈ 0), and let Palioxis
  compute the trapped equilibrium.
- If you need the full mobile+trapped equilibrium simultaneously (e.g. to conserve
  total hydrogen), call `Palioxis.calculate_steady_state!` on each node manually.
"""
function Flopsy.build_equilibrium_ic(
    palioxis_model::Palioxis.MultipleDefectModel,
    adaptor::Flopsy.HotgatesTrappingAdaptor,
    mobile_profile::AbstractVecOrMat,
    T::Real,
)
    n_gas = palioxis_model.n_gas
    n_ne  = palioxis_model.n_ne_species
    nvars = n_gas + n_ne
    nx    = size(adaptor.defects, 2)

    # Normalise to (n_gas, nx) matrix
    mobile_mat = if mobile_profile isa AbstractVector
        n_gas == 1 || throw(ArgumentError(
            "Vector mobile_profile only valid for single-gas models (n_gas=$(n_gas))"
        ))
        reshape(Float64.(mobile_profile), 1, nx)
    else
        Matrix{Float64}(mobile_profile)
    end

    size(mobile_mat) == (n_gas, nx) || throw(DimensionMismatch(
        "mobile_profile must have shape (n_gas=$n_gas, nx=$nx), got $(size(mobile_mat))"
    ))

    u0 = zeros(Float64, nvars * nx)
    U0 = reshape(u0, nvars, nx)

    for ix in 1:nx
        mobile  = mobile_mat[:, ix]
        trapped = Palioxis.set_initial_conditions(palioxis_model, mobile, T)

        for (j, idx) in enumerate(adaptor.mobile_indices)
            U0[idx, ix] = mobile[j]
        end
        for (j, idx) in enumerate(adaptor.trap_indices)
            U0[idx, ix] = trapped[j]
        end
    end

    return u0
end

end # module PalioxisExt
