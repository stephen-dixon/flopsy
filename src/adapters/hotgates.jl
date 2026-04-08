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
"""
struct HotgatesTrappingAdaptor{I,N}
    mobile_indices::Vector{I}
    trap_indices::Vector{I}
    mobile_names::Vector{N}
    trap_names::Vector{N}
    defect_names::Vector{N}
    defects::Matrix{Float64}
end


"""
Reaction operator wrapping a Hotgates-like backend model plus an adaptor.
"""
struct HotgatesReactionOperator{M,A,T} <: AbstractReactionOperator
    model::M
    adaptor::A
    temperature::T
end

supports_rhs(::HotgatesReactionOperator) = true


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
    boundary = nothing,
)
    layout = build_hotgates_variable_layout(adaptor)

    length(diffusion_coefficients) == nvariables(layout) ||
        throw(ArgumentError("diffusion_coefficients must have length $(nvariables(layout))"))

    reaction = HotgatesReactionOperator(model, adaptor, temperature)

    selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)
    diffusion = LinearDiffusionOperator(collect(diffusion_coefficients), selector, nothing)

    return build_rd_model(
        layout = layout,
        mesh = mesh,
        reaction = reaction,
        diffusion = diffusion,
        boundary = boundary,
    )
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
function local_rhs!(du_local, u_local, model, adaptor::HotgatesTrappingAdaptor, ctx, t, ix, T)
    nm = length(adaptor.mobile_indices)
    nt = length(adaptor.trap_indices)
    nd = size(adaptor.defects, 1)

    mobile  = similar(u_local, nm)
    trapped = similar(u_local, nt)
    defects = nd > 0 ? adaptor.defects[:, ix] : similar(u_local, 0)
    dmobile  = similar(u_local, nm)
    dtrapped = similar(u_local, nt)

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


# ------------------------------------------------------------------
# Stub for the real Palioxis backend — implemented in ext/PalioxisExt.jl
# ------------------------------------------------------------------

"""
    build_palioxis_trapping_model(; palioxis_model, mesh, defects, temperature,
                                    diffusion_coefficients)

Build a Flopsy model backed by a real `Palioxis.MultipleDefectModel`.

Requires the `Palioxis` package to be loaded (provided via the `PalioxisExt`
package extension).  Loading `Palioxis` alongside `Flopsy` automatically
activates the extension.

# Arguments
- `palioxis_model` — a constructed `Palioxis.MultipleDefectModel`
- `mesh`           — `Mesh1D` defining the spatial domain
- `defects`        — `Matrix{Float64}` of shape `(n_traps, nx)` with the
                     spatially-varying defect concentration profile.
                     Use `Palioxis.get_defect_concentrations(model, z1, z2)`
                     to build per-node columns.
- `temperature`    — an `AbstractTemperatureProvider`
- `diffusion_coefficients` — length-`n_gas` vector (traps do not diffuse)
"""
function build_palioxis_trapping_model end


"""
    build_equilibrium_ic(palioxis_model, adaptor, mobile_profile, T) -> Vector{Float64}

Build a flat initial-state vector where mobile species concentrations are set from
`mobile_profile` and trapped species are placed at their Palioxis equilibrium values
for temperature `T`.

`mobile_profile` may be:
- A `Vector{Float64}` of length `nx` (single mobile species)
- A `Matrix{Float64}` of shape `(n_gas, nx)` (multiple mobile species)

Calls `Palioxis.set_initial_conditions` per spatial node.
Requires the `Palioxis` package extension to be loaded.

See also `Palioxis.calculate_steady_state!` for finding full equilibrium when
the initial mobile/trapped split is not known a priori.
"""
function build_equilibrium_ic end


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
function hotgates_rates!(dmobile, dtrapped, model::FakeHotgatesModel, mobile, defects, trapped, T)
    fill!(dmobile,  zero(eltype(dmobile)))
    fill!(dtrapped, zero(eltype(dtrapped)))

    length(mobile) == length(trapped) ||
        throw(ArgumentError("FakeHotgatesModel expects equal numbers of mobile and trap variables"))

    @inbounds for j in eachindex(mobile, trapped, dmobile, dtrapped)
        trap_flux   = model.k_trap   * mobile[j] * (1 - trapped[j])
        detrap_flux = model.k_detrap * trapped[j]
        net = trap_flux - detrap_flux

        dmobile[j]  -= net
        dtrapped[j] += net
    end

    return nothing
end
