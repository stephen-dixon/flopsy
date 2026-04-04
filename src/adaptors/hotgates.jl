struct HotgatesReactionOperator{M, A} <: AbstractReactionOperator
    model::M
    adaptor::A
end

supports_rhs(::HotgatesReactionOperator) = true

"""
    rhs!(du, op::HotgatesReactionOperator, u, ctx, t)

Node-local reaction contribution hook for a Hotgates-backed model.

`op.adaptor` is expected to provide:

    local_rhs!(du_local, u_local, model, ctx, t, ix)

where:
- `u_local` is a view of all state variables at node `ix`
- `du_local` is the corresponding output view
"""
function rhs!(du, op::HotgatesReactionOperator, u, ctx::SystemContext, t)
    layout = ctx.aux[:layout]
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
