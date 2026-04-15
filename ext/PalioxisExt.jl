"""
Package extension that wires the real `Palioxis.MultipleDefectModel` into
Flopsy's Hotgates adapter.  Activated automatically when both `Flopsy` and
`Palioxis` are loaded in the same Julia session.
"""
module PalioxisExt

using Flopsy
using Palioxis

# ---------------------------------------------------------------------------
# PalioxisDiffusionCoefficients
# ---------------------------------------------------------------------------

"""
    PalioxisDiffusionCoefficients(model)

`AbstractDiffusionCoefficients` backed by a live `Palioxis.MultipleDefectModel`.

Calls `Palioxis.diffusion_constants(model, T)` at every evaluation, so D is
fully consistent with the Palioxis model parameters and the current temperature.
No pre-evaluation at a reference temperature; no approximation.
"""
struct PalioxisDiffusionCoefficients <: Flopsy.AbstractDiffusionCoefficients
    model::Palioxis.MultipleDefectModel
end

function Flopsy.get_D(c::PalioxisDiffusionCoefficients, ivar::Int, ix::Int, T::Real)
    return Palioxis.diffusion_constants(c.model, T)[ivar]
end

# ---------------------------------------------------------------------------
# hotgates_rates! — dispatch for the real Palioxis backend
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# supports_jacobian — enabled for the Palioxis backend
# ---------------------------------------------------------------------------

function Flopsy.supports_jacobian(op::Flopsy.HotgatesReactionOperator{<:Palioxis.MultipleDefectModel})
    return true
end

# ---------------------------------------------------------------------------
# hotgates_jacobian! — analytic per-node Jacobian via Palioxis C API
# ---------------------------------------------------------------------------

function Flopsy.hotgates_jacobian!(
        J_local::AbstractMatrix,
        model::Palioxis.MultipleDefectModel,
        adaptor::Flopsy.HotgatesTrappingAdaptor,
        mobile::AbstractVector,
        defects::AbstractVector,
        trapped::AbstractVector,
        T::Real
)
    jac = Palioxis.time_derivatives_jacobian(model, defects, mobile, trapped, T)

    n_gas = model.n_gas
    n_ne = model.n_ne_species

    # Unpack blocks — stored Fortran col-major, which matches Julia col-major,
    # so reshape gives the correct matrix directly.
    ddx_dx = reshape(jac.ddx_dx.data, n_gas, n_gas)    # d(mobile_rates)/d(mobile)
    ddx_dyc = reshape(jac.ddx_dyc.data, n_gas, n_ne)     # d(mobile_rates)/d(trapped)
    ddydt_dx = reshape(jac.ddydt_dx.data, n_ne, n_gas)    # d(trapped_rates)/d(mobile)
    ddydt_dyc = reshape(jac.ddydt_dyc.data, n_ne, n_ne)     # d(trapped_rates)/d(trapped)

    m_idxs = adaptor.mobile_indices   # row/col positions in local state
    t_idxs = adaptor.trap_indices

    # Accumulate into J_local (row = output variable, col = input variable)
    @inbounds for (im2, cm) in enumerate(m_idxs), (im1, rm) in enumerate(m_idxs)

        J_local[rm, cm] += ddx_dx[im1, im2]
    end
    @inbounds for (it2, ct) in enumerate(t_idxs), (im1, rm) in enumerate(m_idxs)

        J_local[rm, ct] += ddx_dyc[im1, it2]
    end
    @inbounds for (im2, cm) in enumerate(m_idxs), (it1, rt) in enumerate(t_idxs)

        J_local[rt, cm] += ddydt_dx[it1, im2]
    end
    @inbounds for (it2, ct) in enumerate(t_idxs), (it1, rt) in enumerate(t_idxs)

        J_local[rt, ct] += ddydt_dyc[it1, it2]
    end

    return J_local
end

# ---------------------------------------------------------------------------
# hotgates_rates!
# ---------------------------------------------------------------------------

function Flopsy.hotgates_rates!(
        dmobile::AbstractVector,
        dtrapped::AbstractVector,
        model::Palioxis.MultipleDefectModel,
        mobile::AbstractVector,
        defects::AbstractVector,
        trapped::AbstractVector,
        T::Real
)
    Palioxis.rates!(dmobile, dtrapped, model, mobile, defects, trapped, T)
    return nothing
end

# ---------------------------------------------------------------------------
# build_palioxis_trapping_model
# ---------------------------------------------------------------------------

function Flopsy.build_palioxis_trapping_model(;
        palioxis_model::Palioxis.MultipleDefectModel,
        mesh::Flopsy.Mesh1D,
        defects::AbstractMatrix,
        temperature::Flopsy.AbstractTemperatureProvider,
        left_bc = nothing,
        right_bc = nothing
)
    nx = length(mesh.x)

    size(defects) == (palioxis_model.n_traps, nx) ||
        throw(DimensionMismatch(
            "defects must have shape ($(palioxis_model.n_traps), $nx), " *
            "got $(size(defects))"
        ))

    n_gas = palioxis_model.n_gas
    n_ne = palioxis_model.n_ne_species

    # Build trap_groups from the Palioxis occupancy structure.
    # Each defect type has max_occ[d] fill levels; levels within a type are coupled
    # (multi-occupancy trapping couples adjacent occupancy states).
    trap_occupancies = Palioxis.get_max_trap_occupancy(palioxis_model)
    trap_groups = Vector{Vector{Int}}()
    pos = 1
    for occ in trap_occupancies
        push!(trap_groups, collect(pos:(pos + occ - 1)))
        pos += occ
    end

    adaptor = Flopsy.HotgatesTrappingAdaptor(
        collect(1:n_gas),
        collect((n_gas + 1):(n_gas + n_ne)),
        Palioxis.gas_names(palioxis_model),
        Palioxis.trap_names(palioxis_model),
        Palioxis.defect_names(palioxis_model),
        Matrix{Float64}(defects),
        trap_groups
    )

    # Diffusion coefficients come from Palioxis at each time step — fully T-dependent.
    # The coefficients object wraps the Palioxis model; the selector restricts
    # application to mobile (diffusing) variables only.
    diff_coeffs = PalioxisDiffusionCoefficients(palioxis_model)

    # Flat coefficient vector for LinearDiffusionOperator: the selector only
    # accesses mobile indices so trapped slots are never read, but we need
    # the *diffusion* operator to store the same PalioxisDiffusionCoefficients.
    selector = layout -> Flopsy.variables_with_tag(layout, :diffusion)

    diffusion = Flopsy.LinearDiffusionOperator(diff_coeffs, selector, nothing, temperature)

    boundary = if left_bc !== nothing || right_bc !== nothing
        Flopsy.WeakDirichletBoundaryOperator(selector, diff_coeffs, temperature;
            left = left_bc, right = right_bc)
    else
        nothing
    end

    return Flopsy.build_rd_model(
        layout = Flopsy.build_hotgates_variable_layout(adaptor),
        mesh = mesh,
        reaction = Flopsy.HotgatesReactionOperator(palioxis_model, adaptor, temperature),
        diffusion = diffusion,
        boundary = boundary
    )
end

# ---------------------------------------------------------------------------
# build_equilibrium_ic — mobile profile → Palioxis equilibrium trapped
# ---------------------------------------------------------------------------

function Flopsy.build_equilibrium_ic(
        palioxis_model::Palioxis.MultipleDefectModel,
        model::Flopsy.SystemModel,
        mobile_profile::AbstractVecOrMat,
        T::Real
)
    adaptor = model.operators.reaction.adaptor
    n_gas = palioxis_model.n_gas
    n_ne = palioxis_model.n_ne_species
    nvars = n_gas + n_ne
    nx = size(adaptor.defects, 2)

    mobile_mat = if mobile_profile isa AbstractVector
        n_gas == 1 || throw(ArgumentError(
            "Vector mobile_profile only valid for single-gas models (n_gas=$(n_gas))"
        ))
        reshape(Float64.(mobile_profile), 1, nx)
    else
        Matrix{Float64}(mobile_profile)
    end

    size(mobile_mat) == (n_gas, nx) || throw(DimensionMismatch(
        "mobile_profile must have shape ($n_gas, $nx), got $(size(mobile_mat))"
    ))

    u0 = zeros(Float64, nvars * nx)
    U0 = reshape(u0, nvars, nx)

    for ix in 1:nx
        mobile = mobile_mat[:, ix]
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

# ---------------------------------------------------------------------------
# build_ic_from_total_hydrogen — total H profile → Palioxis steady-state
# ---------------------------------------------------------------------------

function Flopsy.build_ic_from_total_hydrogen(
        palioxis_model::Palioxis.MultipleDefectModel,
        model::Flopsy.SystemModel,
        total_hydrogen::AbstractVector,
        T::Real
)
    adaptor = model.operators.reaction.adaptor
    n_gas = palioxis_model.n_gas
    n_ne = palioxis_model.n_ne_species
    nvars = n_gas + n_ne
    nx = size(adaptor.defects, 2)

    length(total_hydrogen) == nx || throw(DimensionMismatch(
        "total_hydrogen must have length nx=$nx, got $(length(total_hydrogen))"
    ))

    u0 = zeros(Float64, nvars * nx)
    U0 = reshape(u0, nvars, nx)

    for ix in 1:nx
        C_tot = Float64(total_hydrogen[ix])
        defects = adaptor.defects[:, ix]

        # Initial guess: divide evenly across all DOFs, then clamp to valid bounds.
        frac = C_tot / (n_gas + n_ne)
        mobile = fill(frac, n_gas)
        trapped = fill(frac, n_ne)

        Palioxis.ensure_bounds!(palioxis_model, mobile, collect(Float64, defects), trapped)

        # Iterate to the nearest equilibrium point (conserves chemistry, not total H).
        Palioxis.calculate_steady_state!(palioxis_model, mobile, defects, trapped, T)

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
