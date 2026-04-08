struct HotgatesReactionOperator{M,A} <: AbstractReactionOperator
    model::M
    adaptor::A
end

supports_rhs(::HotgatesReactionOperator) = true

"""
Expected adaptor interface:

    local_rhs!(du_local, u_local, model, ctx, t, ix)

where `u_local` and `du_local` are views over all variables at node `ix`.
"""
function rhs!(du, op::HotgatesReactionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx = ctx.nx

    U = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    @inbounds for ix in 1:nx
        u_local = node_view(U, ix)
        du_local = node_view(dU, ix)
        op.adaptor.local_rhs!(du_local, u_local, op.model, ctx, t, ix)
    end

    return du
end
