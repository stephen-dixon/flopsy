struct LinearDiffusionOperator{T,S,B} <: AbstractDiffusionOperator
    coefficients::Vector{T}
    selector::S
    bc::B
end

supports_rhs(::LinearDiffusionOperator) = true

function rhs!(du, op::LinearDiffusionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx = ctx.nx

    U = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    dx = ctx.mesh.dx
    invdx2 = inv(dx * dx)

    vars = op.selector(layout)

    @inbounds for ivar in vars
        D = op.coefficients[ivar]

        dU[ivar, 1] += D * (U[ivar, 2] - U[ivar, 1]) * invdx2

        for ix in 2:(nx - 1)
            dU[ivar, ix] += D * (U[ivar, ix + 1] - 2U[ivar, ix] + U[ivar, ix - 1]) * invdx2
        end

        dU[ivar, nx] += D * (U[ivar, nx - 1] - U[ivar, nx]) * invdx2
    end

    return du
end
