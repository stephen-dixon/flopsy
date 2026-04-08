abstract type AbstractTemperatureProvider end

struct ConstantTemperature{T} <: AbstractTemperatureProvider
    value::T
end

temperature_at(tp::ConstantTemperature, ctx, t, ix) = tp.value


"""
Mapping metadata between Flopsy local node state and a Hotgates-like backend.

Assumes a local solver state vector contains all variables for one spatial node.
"""
struct HotgatesTrappingAdaptor{I,N}
    mobile_indices::Vector{I}
    trap_indices::Vector{I}
    mobile_names::Vector{N}
    trap_names::Vector{N}
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

`diffusion_coefficients` should be length nvariables(layout), with zeros for
non-diffusing variables.
"""
function build_hotgates_trapping_model(;
    mesh::Mesh1D,
    model,
    adaptor::HotgatesTrappingAdaptor,
    temperature::AbstractTemperatureProvider,
    diffusion_coefficients::AbstractVector,
)
    layout = build_hotgates_variable_layout(adaptor)

    length(diffusion_coefficients) == nvariables(layout) ||
        throw(ArgumentError("diffusion_coefficients must have length $(nvariables(layout))"))

    reaction = HotgatesReactionOperator(model, adaptor, temperature)

    selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)
    # diffusion = LinearDiffusionOperator(collect(diffusion_coefficients), selector)
    diffusion = LinearDiffusionOperator(collect(diffusion_coefficients), selector, nothing)

    return build_rd_model(
        layout = layout,
        mesh = mesh,
        reaction = reaction,
        diffusion = diffusion,
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

This is intentionally simple and allocates per node for correctness-first
development. We can replace with reusable caches later.
"""
function local_rhs!(du_local, u_local, model, adaptor::HotgatesTrappingAdaptor, ctx, t, ix, T)
    nm = length(adaptor.mobile_indices)
    nt = length(adaptor.trap_indices)

    mobile = similar(u_local, nm)
    trapped = similar(u_local, nt)
    dmobile = similar(u_local, nm)
    dtrapped = similar(u_local, nt)

    @inbounds for (j, idx) in enumerate(adaptor.mobile_indices)
        mobile[j] = u_local[idx]
    end

    @inbounds for (j, idx) in enumerate(adaptor.trap_indices)
        trapped[j] = u_local[idx]
    end

    hotgates_rates!(dmobile, dtrapped, model, mobile, trapped, T)

    @inbounds for (j, idx) in enumerate(adaptor.mobile_indices)
        du_local[idx] += dmobile[j]
    end

    @inbounds for (j, idx) in enumerate(adaptor.trap_indices)
        du_local[idx] += dtrapped[j]
    end

    return du_local
end


# ------------------------------------------------------------------
# Fake backend for testing the adaptor path before wiring real Hotgates
# ------------------------------------------------------------------

"""
A fake Hotgates-like backend model.

Implements simple trapping/detrapping for any number of mobile and trap variables:
for each pair j,
    mobile_j  <-> trap_j

using:
    trap_flux   = k_trap   * mobile_j * (1 - trap_j)
    detrap_flux = k_detrap * trap_j
"""
struct FakeHotgatesModel{T}
    k_trap::T
    k_detrap::T
end


"""
Hotgates-like in-place rates API expected by the adaptor.

Arguments:
- dmobile, dtrapped: outputs
- model: backend model object
- mobile, trapped: local state vectors
- T: temperature
"""
function hotgates_rates!(dmobile, dtrapped, model::FakeHotgatesModel, mobile, trapped, T)
    fill!(dmobile, zero(eltype(dmobile)))
    fill!(dtrapped, zero(eltype(dtrapped)))

    length(mobile) == length(trapped) ||
        throw(ArgumentError("FakeHotgatesModel expects same number of mobile and trap variables"))

    @inbounds for j in eachindex(mobile, trapped, dmobile, dtrapped)
        trap_flux = model.k_trap * mobile[j] * (1 - trapped[j])
        detrap_flux = model.k_detrap * trapped[j]
        net = trap_flux - detrap_flux

        dmobile[j] -= net
        dtrapped[j] += net
    end

    return nothing
end
