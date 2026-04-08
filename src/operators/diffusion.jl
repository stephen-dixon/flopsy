# ---------------------------------------------------------------------------
# Internal helper — evaluate D for one variable at a node and temperature.
# Dispatch handles both plain vectors (constant) and AbstractDiffusionCoefficients.
# The ix argument is passed through for future spatially-varying D support;
# all current implementations ignore it.
# ---------------------------------------------------------------------------

_eval_D(coeffs::AbstractVector{<:Real}, ivar::Int, ix::Int, T) = Float64(coeffs[ivar])
_eval_D(coeffs::AbstractDiffusionCoefficients, ivar::Int, ix::Int, T) = get_D(coeffs, ivar, ix, T)


# ---------------------------------------------------------------------------
# Helper — variable indices that this operator couples spatially (inter-node).
# Used by build_unsplit_problem to construct the Jacobian sparsity pattern.
# ---------------------------------------------------------------------------

diffusion_variable_indices(::AbstractOperator, layout) = Int[]


# ---------------------------------------------------------------------------
# LinearDiffusionOperator
# ---------------------------------------------------------------------------

"""
    LinearDiffusionOperator(coefficients, selector, bc[, temperature])

Standard second-order finite-difference diffusion operator for 1D domains.

Uses a one-sided stencil at the boundary nodes that corresponds to zero-flux
(Neumann) boundary conditions.  Pair with `DirichletBoundaryOperator` to impose
Dirichlet conditions.

# Arguments
- `coefficients` — diffusion coefficients per variable.  May be:
  - `Vector{<:Real}` — constant, temperature-independent
  - `AbstractDiffusionCoefficients` — e.g. `ArrheniusDiffusion`, `CallableDiffusion`,
    or a library-provided type such as `PalioxisDiffusionCoefficients`
- `selector`  — `f(layout) -> Vector{Int}` returning the variable indices to diffuse.
                Typically `layout -> variables_with_tag(layout, :diffusion)`.
- `bc`        — reserved for future boundary-condition metadata; pass `nothing`.
- `temperature` — `AbstractTemperatureProvider` required when `coefficients` is
                  temperature-dependent; `nothing` for constant coefficients.
"""
struct LinearDiffusionOperator{C, S, B, TP} <: AbstractDiffusionOperator
    coefficients::C
    selector::S
    bc::B
    temperature::TP
end

# Backward-compatible 3-argument constructor (constant D, no temperature provider).
LinearDiffusionOperator(coefficients, selector, bc) =
    LinearDiffusionOperator(coefficients, selector, bc, nothing)

supports_rhs(::LinearDiffusionOperator) = true
supports_jacobian(::LinearDiffusionOperator) = true

diffusion_variable_indices(op::LinearDiffusionOperator, layout) = op.selector(layout)

function rhs!(du, op::LinearDiffusionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    dx     = ctx.mesh.dx
    invdx2 = inv(dx * dx)

    T_val = op.temperature !== nothing ?
        Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN

    U  = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    vars = op.selector(layout)

    @inbounds for ivar in vars
        D = _eval_D(op.coefficients, ivar, 1, T_val)

        dU[ivar, 1] += D * (U[ivar, 2] - U[ivar, 1]) * invdx2

        for ix in 2:(nx - 1)
            dU[ivar, ix] += D * (U[ivar, ix + 1] - 2U[ivar, ix] + U[ivar, ix - 1]) * invdx2
        end

        dU[ivar, nx] += D * (U[ivar, nx - 1] - U[ivar, nx]) * invdx2
    end

    return du
end

function jacobian!(J, op::LinearDiffusionOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    dx     = ctx.mesh.dx
    invdx2 = inv(dx * dx)
    nvars  = nvariables(layout)

    T_val = op.temperature !== nothing ?
        Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN

    vars = op.selector(layout)

    @inbounds for ivar in vars
        D = _eval_D(op.coefficients, ivar, 1, T_val)

        # Node 1 — one-sided Neumann stencil: dU[ivar,1] += D*(U[2]-U[1])*invdx²
        r = ivar
        J[r, r]          += -D * invdx2
        J[r, r + nvars]  +=  D * invdx2

        # Interior nodes
        for ix in 2:(nx - 1)
            r = (ix - 1) * nvars + ivar
            J[r, r - nvars] +=  D * invdx2
            J[r, r]         += -2 * D * invdx2
            J[r, r + nvars] +=  D * invdx2
        end

        # Node nx — one-sided Neumann stencil
        r = (nx - 1) * nvars + ivar
        J[r, r - nvars] +=  D * invdx2
        J[r, r]         += -D * invdx2
    end

    return J
end


# ---------------------------------------------------------------------------
# DirichletBoundaryOperator
# ---------------------------------------------------------------------------

"""
    DirichletBoundaryOperator(selector, coefficients[, temperature]; left=nothing, right=nothing)

Adds ghost-node Dirichlet corrections to boundary nodes for diffusing variables,
on top of the zero-flux (Neumann) stencil in `LinearDiffusionOperator`.

`left` and `right` are callables `f(t) -> concentration_value`, or `nothing` to
leave that boundary at zero-flux.  The correction at boundary node 1 is:

    dU[ivar, 1]  += D(T) * (left(t)  - U[ivar, 1])  / dx²
    dU[ivar, nx] += D(T) * (right(t) - U[ivar, nx]) / dx²

Combined with the Neumann stencil in `LinearDiffusionOperator` this gives the
standard centred Dirichlet stencil `D*(U[2] - 2*U[1] + g)/dx²`.

The `coefficients` argument accepts the same types as `LinearDiffusionOperator`:
a plain `Vector{<:Real}` for constant D, or any `AbstractDiffusionCoefficients`
for temperature-dependent D.  Provide a matching `temperature` provider when
using temperature-dependent coefficients.

# TDS usage — vacuum at both surfaces

```julia
boundary = DirichletBoundaryOperator(
    selector, diffusion_coefficients, temperature;
    left  = t -> 0.0,
    right = t -> 0.0,
)
```

!!! note
    Set the initial condition at boundary nodes to match `left(0)` / `right(0)` to
    avoid a step discontinuity at t = 0.
"""
struct DirichletBoundaryOperator{C, S, TP, L, R} <: AbstractDiffusionOperator
    selector::S
    coefficients::C
    temperature::TP
    left::L    # Union{Nothing, callable f(t) -> value}
    right::R   # Union{Nothing, callable f(t) -> value}
end

function DirichletBoundaryOperator(selector, coefficients, temperature=nothing;
                                    left=nothing, right=nothing)
    return DirichletBoundaryOperator(selector, coefficients, temperature, left, right)
end

supports_rhs(::DirichletBoundaryOperator) = true
supports_jacobian(::DirichletBoundaryOperator) = true

diffusion_variable_indices(op::DirichletBoundaryOperator, layout) = op.selector(layout)

function rhs!(du, op::DirichletBoundaryOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    dx     = ctx.mesh.dx
    invdx2 = inv(dx * dx)

    T_val = op.temperature !== nothing ?
        Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN

    U  = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    vars = op.selector(layout)

    @inbounds for ivar in vars
        D = _eval_D(op.coefficients, ivar, 1, T_val)

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

function jacobian!(J, op::DirichletBoundaryOperator, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    dx     = ctx.mesh.dx
    invdx2 = inv(dx * dx)
    nvars  = nvariables(layout)

    T_val = op.temperature !== nothing ?
        Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN

    vars = op.selector(layout)

    @inbounds for ivar in vars
        D = _eval_D(op.coefficients, ivar, 1, T_val)

        if op.left !== nothing
            r = ivar   # node 1
            J[r, r] += -D * invdx2
        end

        if op.right !== nothing
            r = (nx - 1) * nvars + ivar   # node nx
            J[r, r] += -D * invdx2
        end
    end

    return J
end


# ---------------------------------------------------------------------------
# Surface flux computation — dispatched on diffusion operator type
# ---------------------------------------------------------------------------

"""
    surface_fluxes(op, result) -> Dict{Symbol, NamedTuple}

Compute outward diffusive fluxes at both domain boundaries for each diffusing
variable over all saved times.  Dispatches on the concrete operator type.

Default returns an empty `Dict`.  Implement for concrete `AbstractDiffusionOperator`
subtypes to enable flux reporting.
"""
surface_fluxes(::AbstractDiffusionOperator, result) = Dict{Symbol, NamedTuple}()

function surface_fluxes(op::LinearDiffusionOperator, result::SimulationResult)
    model  = result.model
    layout = model.layout
    mesh   = model.context.mesh
    dx     = mesh.dx
    nx     = model.context.nx
    vars   = op.selector(layout)
    names  = variable_names(layout)
    nt     = length(result.solution.u)
    ctx    = model.context

    out = Dict{Symbol, NamedTuple}()

    for ivar in vars
        left_flux  = zeros(Float64, nt)
        right_flux = zeros(Float64, nt)

        for it in 1:nt
            t     = result.solution.t[it]
            T_val = op.temperature !== nothing ?
                Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN
            D = _eval_D(op.coefficients, ivar, 1, T_val)

            U = state_view(result.solution.u[it], layout, nx)
            left_flux[it]  = D * (U[ivar, 2]      - U[ivar, 1])  / dx
            right_flux[it] = D * (U[ivar, nx - 1] - U[ivar, nx]) / dx
        end

        out[names[ivar]] = (left = left_flux, right = right_flux)
    end

    return out
end
