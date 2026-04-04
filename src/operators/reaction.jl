struct ToyReactionOperator{T} <: AbstractReactionOperator
    rate::T
end

supports_rhs(::ToyReactionOperator) = true

"""
Simple local decay reaction:
    du/dt = -k*u
for every state variable at every node.
"""
function rhs!(du, op::ToyReactionOperator, u, ctx::SystemContext, t)
    @inbounds @. du += -op.rate * u
    return du
end


"""
A simple two-variable local trapping model.

Per node, the local state is assumed to contain:
- c      : a mobile/diffusing concentration-like variable
- theta  : a trap occupancy-like variable

Dynamics:
    dc/dt     += -k_trap * c * (1 - theta) + k_detrap * theta
    dtheta/dt +=  k_trap * c * (1 - theta) - k_detrap * theta

This operator does not assume anything about diffusion.
Diffusion is handled separately by a diffusion operator acting on selected variables.
"""
struct SimpleTrappingReactionOperator{T,I} <: AbstractReactionOperator
    k_trap::T
    k_detrap::T
    mobile_index::I
    trap_index::I
end

supports_rhs(::SimpleTrappingReactionOperator) = true

function rhs!(du, op::SimpleTrappingReactionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx = ctx.nx

    U = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    im = op.mobile_index
    it = op.trap_index

    @inbounds for ix in 1:nx
        c = U[im, ix]
        θ = U[it, ix]

        trap_flux = op.k_trap * c * (1 - θ)
        detrap_flux = op.k_detrap * θ
        net = trap_flux - detrap_flux

        dU[im, ix] -= net
        dU[it, ix] += net
    end

    return du
end
