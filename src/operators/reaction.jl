"""
    ToyReactionOperator(rate)

Simple uniform decay reaction: `du/dt += -rate * u` for every state variable
at every node.  Used for smoke tests and simple validation problems.
"""
struct ToyReactionOperator{T} <: AbstractReactionOperator
    rate::T
end

supports_rhs(::ToyReactionOperator) = true
supports_jacobian(::ToyReactionOperator) = true

function rhs!(du, op::ToyReactionOperator, u, ctx::SystemContext, t)
    @inbounds @. du += -op.rate * u
    return du
end

function jacobian!(J, op::ToyReactionOperator, u, ctx::SystemContext, t)
    n = length(u)
    @inbounds for i in 1:n
        J[i, i] += -op.rate
    end
    return J
end

jacobian_node_sparsity(::ToyReactionOperator, layout) =
    Set{Tuple{Int,Int}}((i, i) for i in 1:nvariables(layout))


"""
    SimpleTrappingReactionOperator(k_trap, k_detrap, mobile_index, trap_index)

Two-variable local trapping model.

Per node, the local state contains a mobile concentration `c` and a trap
occupancy `θ` (dimensionless, 0–1).  Dynamics:

    dc/dt     += -k_trap * c * (1 - θ) + k_detrap * θ
    dθ/dt     +=  k_trap * c * (1 - θ) - k_detrap * θ

Diffusion is handled separately by a `LinearDiffusionOperator` on the mobile
variable.  `mobile_index` and `trap_index` are 1-based variable indices in the
layout.
"""
struct SimpleTrappingReactionOperator{T,I} <: AbstractReactionOperator
    k_trap::T
    k_detrap::T
    mobile_index::I
    trap_index::I
end

supports_rhs(::SimpleTrappingReactionOperator) = true
supports_jacobian(::SimpleTrappingReactionOperator) = true

function rhs!(du, op::SimpleTrappingReactionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx = ctx.nx

    U  = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    im = op.mobile_index
    it = op.trap_index

    @inbounds for ix in 1:nx
        c = U[im, ix]
        θ = U[it, ix]

        trap_flux   = op.k_trap   * c * (1 - θ)
        detrap_flux = op.k_detrap * θ
        net = trap_flux - detrap_flux

        dU[im, ix] -= net
        dU[it, ix] += net
    end

    return du
end

function jacobian!(J, op::SimpleTrappingReactionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    nvars  = nvariables(layout)

    U  = state_view(u, layout, nx)
    im = op.mobile_index
    it = op.trap_index

    @inbounds for ix in 1:nx
        c = U[im, ix]
        θ = U[it, ix]

        row_m = (ix - 1) * nvars + im
        row_t = (ix - 1) * nvars + it

        # Partial derivatives of net = k_trap*c*(1-θ) - k_detrap*θ
        # dc/dt += -net  →  d(dc/dt)/d(c) = -k_trap*(1-θ), d(dc/dt)/d(θ) = k_trap*c + k_detrap
        # dθ/dt += +net  →  d(dθ/dt)/d(c) = k_trap*(1-θ),  d(dθ/dt)/d(θ) = -k_trap*c - k_detrap

        J[row_m, row_m] += -op.k_trap * (1 - θ)
        J[row_m, row_t] +=  op.k_trap * c + op.k_detrap

        J[row_t, row_m] +=  op.k_trap * (1 - θ)
        J[row_t, row_t] += -op.k_trap * c - op.k_detrap
    end

    return J
end

function jacobian_node_sparsity(op::SimpleTrappingReactionOperator, layout)
    im = op.mobile_index
    it = op.trap_index
    return Set{Tuple{Int,Int}}([(im, im), (im, it), (it, im), (it, it)])
end
