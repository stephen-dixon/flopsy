"""
Mapping metadata between Flopsy's local node state and a Hotgates-like backend.

Fields:
- `mobile_indices`  — indices of mobile species in the local Flopsy state vector
- `trap_indices`    — indices of trapped species in the local Flopsy state vector
- `mobile_names`    — species names for mobile variables (metadata / postprocessing)
- `trap_names`      — species names for trapped variables
- `defect_names`    — species names for defect sites (not in state vector)
- `defects`         — spatially-varying defect concentration profile,
                      shape `(n_defects, nx)`.  Column `ix` is the defect vector
                      at node `ix`.  Use `zeros(Float64, 0, nx)` for backends
                      that do not require defect data (e.g. FakeHotgatesModel).
- `trap_groups`     — groups of positions (1-based) within `trap_indices` that belong to
                      the same defect type.  Used to determine the per-node Jacobian sparsity
                      pattern: trap levels within a group are coupled (tridiagonal block);
                      levels in different groups are not coupled.
                      Default (6-arg constructor): each trap in its own singleton group.
"""
struct HotgatesTrappingAdaptor{I, N}
    mobile_indices::Vector{I}
    trap_indices::Vector{I}
    mobile_names::Vector{N}
    trap_names::Vector{N}
    defect_names::Vector{N}
    defects::Matrix{Float64}
    trap_groups::Vector{Vector{Int}}
end

"""
    HotgatesTrappingAdaptor(mobile, trap, mob_names, trap_names, def_names, defects)

Backward-compatible 6-argument constructor.  Sets `trap_groups` so each trap species
is its own singleton group (no intra-group coupling assumed — fully diagonal trap-trap
Jacobian block).  Override by supplying the 7-argument form when defect types have
multiple fill levels that are coupled (e.g. Palioxis multi-occupancy trapping).
"""
function HotgatesTrappingAdaptor(mobile_indices, trap_indices, mobile_names,
        trap_names, defect_names, defects)
    nt = length(trap_indices)
    trap_groups = [[i] for i in 1:nt]
    return HotgatesTrappingAdaptor(mobile_indices, trap_indices, mobile_names,
        trap_names, defect_names, defects, trap_groups)
end

"""
Reaction operator wrapping a Hotgates-like backend model plus an adaptor.
"""
struct HotgatesReactionOperator{M, A, T} <: AbstractReactionOperator
    model::M
    adaptor::A
    temperature::T
end

struct HotgatesWorkspace{T}
    mobile::Vector{T}
    trapped::Vector{T}
    dmobile::Vector{T}
    dtrapped::Vector{T}
    mobile_fd::Vector{T}
    trapped_fd::Vector{T}
    dmobile_fd::Vector{T}
    dtrapped_fd::Vector{T}
    J_local::Matrix{T}
end

function HotgatesWorkspace(nmobile::Int, ntrapped::Int, nvars::Int)
    return HotgatesWorkspace(
        zeros(Float64, nmobile),
        zeros(Float64, ntrapped),
        zeros(Float64, nmobile),
        zeros(Float64, ntrapped),
        zeros(Float64, nmobile),
        zeros(Float64, ntrapped),
        zeros(Float64, nmobile),
        zeros(Float64, ntrapped),
        zeros(Float64, nvars, nvars)
    )
end

supports_rhs(::HotgatesReactionOperator) = true

# jacobian! support is backend-dependent.  The Palioxis extension overrides
# supports_jacobian for HotgatesReactionOperator{<:Palioxis.MultipleDefectModel}
# and implements hotgates_jacobian! via time_derivatives_jacobian.
supports_jacobian(::HotgatesReactionOperator) = false

"""
    jacobian_node_sparsity(op::HotgatesReactionOperator, layout)

Per-node Jacobian sparsity derived from the adaptor's `trap_groups`:
- Mobile-mobile block: fully dense
- Mobile-trap and trap-mobile blocks: fully dense
- Trap-trap block: tridiagonal within each group (adjacent occupancy levels couple)

`trap_groups` is set at adaptor construction time — use the 7-argument
`HotgatesTrappingAdaptor` constructor with defect-type groupings for Palioxis
multi-occupancy models, or leave the default (each trap its own singleton group)
for backends with independent trap species.
"""
function jacobian_node_sparsity(op::HotgatesReactionOperator, layout)
    adaptor = op.adaptor
    m_idxs = adaptor.mobile_indices
    t_idxs = adaptor.trap_indices
    entries = Set{Tuple{Int, Int}}()

    # Mobile-mobile block (dense)
    for cm in m_idxs, rm in m_idxs

        push!(entries, (rm, cm))
    end

    # Mobile←trap and trap←mobile blocks (dense)
    for ct in t_idxs, rm in m_idxs

        push!(entries, (rm, ct))
    end
    for cm in m_idxs, rt in t_idxs

        push!(entries, (rt, cm))
    end

    # Trap-trap: tridiagonal within each group (adjacent levels are coupled)
    for group in adaptor.trap_groups
        gvars = t_idxs[group]  # actual variable indices for this group
        for k in eachindex(gvars)
            push!(entries, (gvars[k], gvars[k]))            # diagonal
            if k > 1
                push!(entries, (gvars[k], gvars[k - 1])) # sub-diagonal
                push!(entries, (gvars[k - 1], gvars[k]))   # super-diagonal
            end
        end
    end

    return entries
end

"""
Construct a Flopsy VariableLayout from adaptor metadata.

Mobile variables are tagged with :reaction and :diffusion.
Trap variables are tagged with :reaction only.
"""
function build_hotgates_variable_layout(adaptor::HotgatesTrappingAdaptor)
    vars = VariableInfo[]

    for name in adaptor.mobile_names
        push!(vars, VariableInfo(Symbol(name), :mobile, Set([:reaction, :diffusion])))
    end

    for name in adaptor.trap_names
        push!(vars, VariableInfo(Symbol(name), :trap, Set([:reaction])))
    end

    return VariableLayout(vars)
end

"""
Build a generic reaction-diffusion model using a Hotgates-like reaction backend.

`diffusion_coefficients` must have length `nvariables(layout)`, with zeros for
non-diffusing variables (typically all trap species).
"""
function build_hotgates_trapping_model(;
        mesh::Mesh1D,
        model,
        adaptor::HotgatesTrappingAdaptor,
        temperature::AbstractTemperatureProvider,
        diffusion_coefficients::AbstractVector,
        boundary = nothing
)
    layout = build_hotgates_variable_layout(adaptor)

    length(diffusion_coefficients) == nvariables(layout) ||
        throw(ArgumentError("diffusion_coefficients must have length $(nvariables(layout))"))

    reaction = HotgatesReactionOperator(model, adaptor, temperature)

    selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)
    diffusion = LinearDiffusionOperator(collect(diffusion_coefficients), selector, nothing)
    scratch = (
        hotgates = HotgatesWorkspace(
        length(adaptor.mobile_indices), length(adaptor.trap_indices), nvariables(layout)),
    )

    return build_rd_model(
        layout = layout,
        mesh = mesh,
        reaction = reaction,
        diffusion = diffusion,
        boundary = boundary,
        scratch = scratch
    )
end

"""
    hotgates_jacobian!(J_local, model, adaptor, mobile, defects, trapped, T)

In-place per-node reaction Jacobian API.  `J_local` is the `(nvars × nvars)`
local Jacobian block for the current node (row = output variable, col = input
variable).  The default (FakeHotgatesModel) falls back to finite differences
by perturbing each input.  The Palioxis extension overrides this for
`Palioxis.MultipleDefectModel` using the analytic `time_derivatives_jacobian`.
"""
function hotgates_jacobian!(J_local, model, adaptor::HotgatesTrappingAdaptor,
        mobile, defects, trapped, T, workspace::HotgatesWorkspace)
    # Default: finite-difference fallback using hotgates_rates!
    nm = length(adaptor.mobile_indices)
    nt = length(adaptor.trap_indices)
    fill!(J_local, 0.0)
    fill!(workspace.mobile_fd, 0.0)
    fill!(workspace.trapped_fd, 0.0)
    fill!(workspace.dmobile_fd, 0.0)
    fill!(workspace.dtrapped_fd, 0.0)

    hotgates_rates!(
        workspace.dmobile_fd, workspace.dtrapped_fd, model, mobile, defects, trapped, T)

    eps = 1e-7
    copyto!(workspace.mobile_fd, mobile)
    copyto!(workspace.trapped_fd, trapped)

    for (j, idx) in enumerate(adaptor.mobile_indices)
        h = eps * max(1.0, abs(mobile[j]))
        workspace.mobile_fd[j] += h
        fill!(workspace.dmobile, 0.0)
        fill!(workspace.dtrapped, 0.0)
        hotgates_rates!(workspace.dmobile, workspace.dtrapped, model,
            workspace.mobile_fd, defects, workspace.trapped_fd, T)
        for i in eachindex(adaptor.mobile_indices)
            J_local[adaptor.mobile_indices[i], idx] += (workspace.dmobile[i] -
                                                        workspace.dmobile_fd[i]) / h
        end
        for i in eachindex(adaptor.trap_indices)
            J_local[adaptor.trap_indices[i], idx] += (workspace.dtrapped[i] -
                                                      workspace.dtrapped_fd[i]) / h
        end
        workspace.mobile_fd[j] -= h
    end

    for (j, idx) in enumerate(adaptor.trap_indices)
        h = eps * max(1.0, abs(trapped[j]))
        workspace.trapped_fd[j] += h
        fill!(workspace.dmobile, 0.0)
        fill!(workspace.dtrapped, 0.0)
        hotgates_rates!(workspace.dmobile, workspace.dtrapped, model,
            workspace.mobile_fd, defects, workspace.trapped_fd, T)
        for i in eachindex(adaptor.mobile_indices)
            J_local[adaptor.mobile_indices[i], idx] += (workspace.dmobile[i] -
                                                        workspace.dmobile_fd[i]) / h
        end
        for i in eachindex(adaptor.trap_indices)
            J_local[adaptor.trap_indices[i], idx] += (workspace.dtrapped[i] -
                                                      workspace.dtrapped_fd[i]) / h
        end
        workspace.trapped_fd[j] -= h
    end

    return J_local
end

"""
Main Flopsy RHS bridge for Hotgates-like backends.
"""
function rhs!(du, op::HotgatesReactionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx = ctx.nx

    U = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    @inbounds for ix in 1:nx
        u_local = node_view(U, ix)
        du_local = node_view(dU, ix)
        T = temperature_at(op.temperature, ctx, t, ix)
        local_rhs!(du_local, u_local, op.model, op.adaptor, ctx, t, ix, T)
    end

    return du
end

"""
Default local bridge for trapping-style Hotgates adaptors.

Extracts mobile and trapped concentrations from the node state, calls
`hotgates_rates!` on the backend model, and accumulates the result back into
the node residual.  Defect concentrations are taken from `adaptor.defects[:, ix]`.
"""
function local_rhs!(
        du_local, u_local, model, adaptor::HotgatesTrappingAdaptor, ctx, t, ix, T)
    nd = size(adaptor.defects, 1)
    workspace = _hotgates_workspace(ctx.scratch, adaptor, length(u_local))
    mobile = workspace.mobile
    trapped = workspace.trapped
    dmobile = workspace.dmobile
    dtrapped = workspace.dtrapped
    defects = nd > 0 ? adaptor.defects[:, ix] : view(adaptor.defects, :, ix)

    fill!(dmobile, 0.0)
    fill!(dtrapped, 0.0)

    @inbounds for (j, idx) in enumerate(adaptor.mobile_indices)
        mobile[j] = u_local[idx]
    end

    @inbounds for (j, idx) in enumerate(adaptor.trap_indices)
        trapped[j] = u_local[idx]
    end

    hotgates_rates!(dmobile, dtrapped, model, mobile, defects, trapped, T)

    @inbounds for (j, idx) in enumerate(adaptor.mobile_indices)
        du_local[idx] += dmobile[j]
    end

    @inbounds for (j, idx) in enumerate(adaptor.trap_indices)
        du_local[idx] += dtrapped[j]
    end

    return du_local
end

"""
Analytic Jacobian bridge for Hotgates-like backends.

Fills the block-diagonal reaction Jacobian in `J`.  Each node contributes an
`nvars × nvars` dense block.  Only called when `supports_jacobian` returns `true`
for the operator (requires the backend to support analytic Jacobians, e.g. Palioxis).
"""
function jacobian!(J, op::HotgatesReactionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx = ctx.nx
    nvars = nvariables(layout)

    U = state_view(u, layout, nx)

    @inbounds for ix in 1:nx
        u_local = node_view(U, ix)
        T = temperature_at(op.temperature, ctx, t, ix)
        nd = size(op.adaptor.defects, 1)
        workspace = _hotgates_workspace(ctx.scratch, op.adaptor, nvars)
        mobile = workspace.mobile
        trapped = workspace.trapped
        defects = nd > 0 ? op.adaptor.defects[:, ix] : Float64[]

        offset = (ix - 1) * nvars
        for (j, idx) in enumerate(op.adaptor.mobile_indices)
            mobile[j] = u_local[idx]
        end
        for (j, idx) in enumerate(op.adaptor.trap_indices)
            trapped[j] = u_local[idx]
        end
        J_local = workspace.J_local
        hotgates_jacobian!(
            J_local, op.model, op.adaptor, mobile, defects, trapped, T, workspace)

        for c in 1:nvars, r in 1:nvars

            J[offset + r, offset + c] += J_local[r, c]
        end
    end

    return J
end

function _hotgates_workspace(scratch::AbstractDict, adaptor::HotgatesTrappingAdaptor, nvars::Int)
    return get!(scratch, :hotgates) do
        HotgatesWorkspace(length(adaptor.mobile_indices), length(adaptor.trap_indices), nvars)
    end
end

function _hotgates_workspace(scratch::NamedTuple, adaptor::HotgatesTrappingAdaptor, nvars::Int)
    if hasproperty(scratch, :hotgates)
        return getproperty(scratch, :hotgates)
    end
    return HotgatesWorkspace(length(adaptor.mobile_indices), length(adaptor.trap_indices), nvars)
end

# ------------------------------------------------------------------
# Stub for the real Palioxis backend — implemented in ext/PalioxisExt.jl
# ------------------------------------------------------------------

"""
    build_palioxis_trapping_model(; palioxis_model, mesh, defects, temperature,
                                    left_bc=nothing, right_bc=nothing)

Build a Flopsy model backed by a real `Palioxis.MultipleDefectModel`.

Requires the `Palioxis` package extension (`using Palioxis` alongside `using Flopsy`).
Diffusion coefficients are queried from Palioxis at every time step via
`Palioxis.diffusion_constants(model, T)` — fully temperature-dependent D at no
extra implementation cost.

# Arguments
- `palioxis_model` — a constructed `Palioxis.MultipleDefectModel`
- `mesh`           — `Mesh1D` defining the spatial domain
- `defects`        — `Matrix{Float64}` of shape `(n_traps, nx)`. Column `ix` is the
                     defect concentration vector at node `ix`.
- `temperature`    — an `AbstractTemperatureProvider`
- `left_bc`        — callable `f(t) -> value` for left surface concentration, or
                     `nothing` (zero-flux / Neumann)
- `right_bc`       — callable `f(t) -> value` for right surface concentration, or
                     `nothing` (zero-flux / Neumann)
"""
function build_palioxis_trapping_model end

"""
    build_equilibrium_ic(palioxis_model, model, mobile_profile, T) -> Vector{Float64}

Build an initial-state vector with mobile species set from `mobile_profile` and
trapped species at their Palioxis equilibrium values for temperature `T`.

`mobile_profile` may be a `Vector` of length `nx` (single mobile species) or a
`Matrix` of shape `(n_gas, nx)`.  Calls `Palioxis.set_initial_conditions` per node.

Requires the `Palioxis` package extension.
"""
function build_equilibrium_ic end

"""
    build_ic_from_total_hydrogen(palioxis_model, model, total_hydrogen, T) -> Vector{Float64}

Build an initial-state vector from a total hydrogen concentration profile by finding
the equilibrium partition between mobile and trapped species at temperature `T`.

`total_hydrogen` is a `Vector` of length `nx`.  Per node:
1. Distribute `total_hydrogen[ix]` as an initial guess across mobile and trapped DOFs.
2. Clamp to valid bounds with `Palioxis.ensure_bounds!`.
3. Call `Palioxis.calculate_steady_state!` to converge to the nearest equilibrium.

!!! note
    `calculate_steady_state!` does not enforce strict conservation of total hydrogen.
    For typical TDS conditions (strong trapping at low T) the result is physically
    accurate; verify with a mass-balance check if precision is critical.

Requires the `Palioxis` package extension.
"""
function build_ic_from_total_hydrogen end

# ------------------------------------------------------------------
# Fake backend for testing the adaptor path before wiring real Palioxis
# ------------------------------------------------------------------

"""
A fake Hotgates-like backend model for testing without the real Palioxis library.

Implements simple trapping/detrapping for any number of mobile/trap pairs:
    trap_flux   = k_trap   * mobile_j * (1 - trap_j)
    detrap_flux = k_detrap * trapped_j
    net         = trap_flux - detrap_flux
"""
struct FakeHotgatesModel{T}
    k_trap::T
    k_detrap::T
end

"""
    hotgates_rates!(dmobile, dtrapped, model, mobile, defects, trapped, T)

In-place reaction rates API used by `local_rhs!`.

`dmobile` and `dtrapped` are output buffers (written, not accumulated).
`defects` is the defect concentration vector for the current node (may be
length-0 for backends that do not use it).
"""
function hotgates_rates!(
        dmobile, dtrapped, model::FakeHotgatesModel, mobile, defects, trapped, T)
    fill!(dmobile, zero(eltype(dmobile)))
    fill!(dtrapped, zero(eltype(dtrapped)))

    length(mobile) == length(trapped) ||
        throw(ArgumentError("FakeHotgatesModel expects equal numbers of mobile and trap variables"))

    @inbounds for j in eachindex(mobile, trapped, dmobile, dtrapped)
        trap_flux = model.k_trap * mobile[j] * (1 - trapped[j])
        detrap_flux = model.k_detrap * trapped[j]
        net = trap_flux - detrap_flux

        dmobile[j] -= net
        dtrapped[j] += net
    end

    return nothing
end
