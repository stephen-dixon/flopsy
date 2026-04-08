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


"""
    DirichletBoundaryOperator(selector, coefficients; left=nothing, right=nothing)

Applies time-varying Dirichlet boundary conditions to diffusing variables by adding
ghost-node corrections on top of the zero-flux (Neumann) stencil in
`LinearDiffusionOperator`.

`left` and `right` are callables `f(t) -> concentration_value`, or `nothing` to
leave that boundary at zero-flux.  The correction at boundary node 1 is:

    dU[ivar, 1]  += D * (left(t)  - U[ivar, 1])  / dx²
    dU[ivar, nx] += D * (right(t) - U[ivar, nx]) / dx²

Combined with the Neumann stencil already in `LinearDiffusionOperator` this gives
the standard centred Dirichlet stencil `D*(U[2] - 2*U[1] + g)/dx²`.

# TDS usage

For a desorption experiment (both surfaces held at vacuum, i.e. zero concentration):

```julia
boundary = DirichletBoundaryOperator(
    selector,
    diffusion_coefficients;
    left  = t -> 0.0,
    right = t -> 0.0,
)
```

For a combined implantation+desorption run with a time-varying surface source:

```julia
boundary = DirichletBoundaryOperator(
    selector,
    diffusion_coefficients;
    left  = t -> t < t_implant ? source_conc(t) : 0.0,
    right = t -> 0.0,
)
```

!!! note
    Set the initial condition at boundary nodes to match `left(0)` / `right(0)` to
    avoid a discontinuity at t=0.
"""
struct DirichletBoundaryOperator{S,C,L,R} <: AbstractDiffusionOperator
    selector::S
    coefficients::C
    left::L   # Union{Nothing, callable f(t)->value}
    right::R  # Union{Nothing, callable f(t)->value}
end

function DirichletBoundaryOperator(selector, coefficients; left=nothing, right=nothing)
    return DirichletBoundaryOperator(selector, coefficients, left, right)
end

supports_rhs(::DirichletBoundaryOperator) = true

function rhs!(du, op::DirichletBoundaryOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx = ctx.nx

    U = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    dx = ctx.mesh.dx
    invdx2 = inv(dx * dx)

    vars = op.selector(layout)

    @inbounds for ivar in vars
        D = op.coefficients[ivar]

        if op.left !== nothing
            g = op.left(t)
            dU[ivar, 1] += D * (g - U[ivar, 1]) * invdx2
        end

        if op.right !== nothing
            g = op.right(t)
            dU[ivar, nx] += D * (g - U[ivar, nx]) * invdx2
        end
    end

    return du
end
