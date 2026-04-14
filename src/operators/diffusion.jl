# ---------------------------------------------------------------------------
# Internal helper — evaluate D for one variable at a node and temperature.
# ---------------------------------------------------------------------------

_eval_D(coeffs::AbstractVector{<:Real}, ivar::Int, ix::Int, T) = Float64(coeffs[ivar])
_eval_D(coeffs::AbstractDiffusionCoefficients, ivar::Int, ix::Int, T) = get_D(coeffs, ivar, ix, T)


# ---------------------------------------------------------------------------
# Helper — variable indices that this operator couples spatially (inter-node).
# Used by build_unsplit_problem to construct the Jacobian sparsity pattern.
# ---------------------------------------------------------------------------

diffusion_variable_indices(::AbstractOperator, layout) = Int[]


# ===========================================================================
# LinearDiffusionOperator
# ===========================================================================

"""
    LinearDiffusionOperator(coefficients, selector, bc[, temperature])

Standard second-order finite-difference diffusion operator for 1D domains.

Uses a one-sided stencil at the boundary nodes that corresponds to zero-flux
(Neumann) boundary conditions.  Pair with a Dirichlet boundary operator to
impose Dirichlet conditions.

# Arguments
- `coefficients` — diffusion coefficients per variable.
- `selector`  — `f(layout) -> Vector{Int}` returning diffusing variable indices.
- `bc`        — reserved; pass `nothing`.
- `temperature` — `AbstractTemperatureProvider` for temperature-dependent D.
"""
struct LinearDiffusionOperator{C, S, B, TP} <: AbstractDiffusionOperator
    coefficients::C
    selector::S
    bc::B
    temperature::TP
end

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

        r = ivar
        J[r, r]          += -D * invdx2
        J[r, r + nvars]  +=  D * invdx2

        for ix in 2:(nx - 1)
            r = (ix - 1) * nvars + ivar
            J[r, r - nvars] +=  D * invdx2
            J[r, r]         += -2 * D * invdx2
            J[r, r + nvars] +=  D * invdx2
        end

        r = (nx - 1) * nvars + ivar
        J[r, r - nvars] +=  D * invdx2
        J[r, r]         += -D * invdx2
    end

    return J
end


# ===========================================================================
# WeakDirichletBoundaryOperator  (ghost-node correction, formerly DirichletBoundaryOperator)
# ===========================================================================

"""
    WeakDirichletBoundaryOperator(selector, coefficients[, temperature]; left=nothing, right=nothing)

Ghost-node weak Dirichlet boundary condition.  Adds a correction to boundary
nodes on top of the zero-flux (Neumann) stencil in `LinearDiffusionOperator`:

    dU[ivar, 1]  += D(T) * (left(t)  - U[ivar, 1])  / dx²
    dU[ivar, nx] += D(T) * (right(t) - U[ivar, nx]) / dx²

Combined with the Neumann stencil this gives the standard centred Dirichlet
stencil `D*(U[2] - 2*U[1] + g)/dx²`.  Both sides are controlled by a single
operator; use `left=nothing` or `right=nothing` to leave a surface at zero-flux.

For per-side independent BC type control, use `DirichletBoundaryOperator` instead,
which is parameterised by method and applies to one side at a time.
"""
struct WeakDirichletBoundaryOperator{C, S, TP, L, R} <: AbstractDiffusionOperator
    selector::S
    coefficients::C
    temperature::TP
    left::L    # Union{Nothing, callable f(t) -> value}
    right::R   # Union{Nothing, callable f(t) -> value}
end

function WeakDirichletBoundaryOperator(selector, coefficients, temperature=nothing;
                                        left=nothing, right=nothing)
    return WeakDirichletBoundaryOperator(selector, coefficients, temperature, left, right)
end

supports_rhs(::WeakDirichletBoundaryOperator) = true
supports_jacobian(::WeakDirichletBoundaryOperator) = true

diffusion_variable_indices(op::WeakDirichletBoundaryOperator, layout) = op.selector(layout)

function rhs!(du, op::WeakDirichletBoundaryOperator, u, ctx::SystemContext, t)
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

function jacobian!(J, op::WeakDirichletBoundaryOperator, u, ctx::SystemContext, t)
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
            r = ivar
            J[r, r] += -D * invdx2
        end

        if op.right !== nothing
            r = (nx - 1) * nvars + ivar
            J[r, r] += -D * invdx2
        end
    end

    return J
end


# ===========================================================================
# Strong DirichletBoundaryOperator — four enforcement methods, per side
# ===========================================================================

"""
    PenaltyMethod(lambda=1e10)

Strong Dirichlet via penalty/relaxation.  Adds `λ*(g(t) - u[boundary])` to the
boundary node's derivative.  Large λ forces fast relaxation to the target value.
Compatible with all ODE solvers but adds stiffness proportional to λ.
"""
struct PenaltyMethod
    lambda::Float64
    PenaltyMethod(λ=1e10) = new(Float64(λ))
end

"""
    MassMatrixMethod()

Strong Dirichlet via a singular mass matrix (DAE formulation).
The boundary node's equation becomes the algebraic constraint `u[boundary] = g(t)`.
Requires a DAE-capable solver such as `Rodas5P()`.
"""
struct MassMatrixMethod end

"""
    CallbackMethod()

Strong Dirichlet via a `DiscreteCallback` that resets the boundary node to
`g(t)` after every accepted step.  The `rhs!` contribution is zero (no-op).
Compatible with any solver; slightly less accurate than penalty/mass-matrix
within a step but exact at step boundaries.
"""
struct CallbackMethod end

"""
    EliminatedMethod()

Strong Dirichlet by eliminating the boundary DOF from the ODE dynamics.
Cancels the Neumann stencil contribution added by `LinearDiffusionOperator` at
the boundary node, replaces it with `g'(t)` (computed by finite difference),
and corrects the interior-neighbor stencil so it sees the exact BC value.

Recommended for constant or slowly varying BCs (e.g. vacuum: g=0).
For rapidly time-varying BCs the small finite-difference error in g'(t) is
negligible compared to time-stepping error.
"""
struct EliminatedMethod end

"""
    DirichletBoundaryOperator(side, bc_fn, selector[, coefficients, temperature]; method=PenaltyMethod())

Strong Dirichlet boundary condition applied to one domain side.

# Arguments
- `side`        — `:left` (node 1) or `:right` (node nx)
- `bc_fn`       — callable `f(t) -> concentration_value`
- `selector`    — `f(layout) -> Vector{Int}` selecting diffusing variable indices
- `coefficients` — diffusion coefficients (required for `EliminatedMethod`; ignored otherwise)
- `temperature` — temperature provider (required for `EliminatedMethod` with temperature-dependent D)

# Keyword arguments
- `method` — one of:
  - `PenaltyMethod(λ=1e10)` *(default)* — stiff relaxation `du[boundary] += λ*(g-u)`
  - `MassMatrixMethod()` — DAE constraint; requires `Rodas5P()` or similar
  - `CallbackMethod()` — discrete callback pins boundary after each step
  - `EliminatedMethod()` — cancels Neumann stencil, sets `du[boundary] = g'(t)`

# Example — independent BC types on left and right
```julia
selector = layout -> variables_with_tag(layout, :diffusion)

# Left: strong penalty, Right: weak ghost-node
left_bc  = DirichletBoundaryOperator(:left,  t -> C_surface, selector; method=PenaltyMethod(1e12))
right_bc = WeakDirichletBoundaryOperator(selector, coeffs, temp; right = t -> 0.0)

# Compose into model boundary slot via OperatorSum
boundary = OperatorSum((left_bc, right_bc))
```
"""
struct DirichletBoundaryOperator{M, S, C, TP, G} <: AbstractDiffusionOperator
    method::M
    side::Symbol      # :left or :right
    bc_fn::G
    selector::S
    coefficients::C   # needed for EliminatedMethod; nothing otherwise
    temperature::TP
end

function DirichletBoundaryOperator(side::Symbol, bc_fn, selector,
                                    coefficients=nothing, temperature=nothing;
                                    method=PenaltyMethod())
    side in (:left, :right) || throw(ArgumentError("side must be :left or :right, got :$side"))
    return DirichletBoundaryOperator(method, side, bc_fn, selector, coefficients, temperature)
end

diffusion_variable_indices(op::DirichletBoundaryOperator, layout) = op.selector(layout)


# --- supports_* traits per method ---

supports_rhs(::DirichletBoundaryOperator) = true
supports_jacobian(::DirichletBoundaryOperator{<:PenaltyMethod})    = true
supports_jacobian(::DirichletBoundaryOperator{<:MassMatrixMethod}) = true
supports_jacobian(::DirichletBoundaryOperator{<:CallbackMethod})   = true
supports_jacobian(::DirichletBoundaryOperator{<:EliminatedMethod}) = true

supports_mass_matrix(op::DirichletBoundaryOperator{<:MassMatrixMethod}) = true


# ---------------------------------------------------------------------------
# PenaltyMethod  rhs! / jacobian!
# ---------------------------------------------------------------------------

function rhs!(du, op::DirichletBoundaryOperator{<:PenaltyMethod}, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    λ      = op.method.lambda

    U  = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    vars = op.selector(layout)
    ix   = op.side === :left ? 1 : nx

    @inbounds for ivar in vars
        g = op.bc_fn(t)
        dU[ivar, ix] += λ * (g - U[ivar, ix])
    end

    return du
end

function jacobian!(J, op::DirichletBoundaryOperator{<:PenaltyMethod}, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    nvars  = nvariables(layout)
    λ      = op.method.lambda

    vars = op.selector(layout)
    ix   = op.side === :left ? 1 : nx

    @inbounds for ivar in vars
        r = (ix - 1) * nvars + ivar
        J[r, r] += -λ
    end

    return J
end


# ---------------------------------------------------------------------------
# MassMatrixMethod  rhs! / jacobian! / mass_matrix
# ---------------------------------------------------------------------------

function rhs!(du, op::DirichletBoundaryOperator{<:MassMatrixMethod}, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx

    U  = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    vars = op.selector(layout)
    ix   = op.side === :left ? 1 : nx

    # Algebraic residual: f[boundary] = g(t) - u[boundary]  (→ u[boundary] = g(t))
    @inbounds for ivar in vars
        g = op.bc_fn(t)
        dU[ivar, ix] += g - U[ivar, ix]
    end

    return du
end

function jacobian!(J, op::DirichletBoundaryOperator{<:MassMatrixMethod}, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    nvars  = nvariables(layout)

    vars = op.selector(layout)
    ix   = op.side === :left ? 1 : nx

    @inbounds for ivar in vars
        r = (ix - 1) * nvars + ivar
        J[r, r] += -1.0
    end

    return J
end

"""
    mass_matrix(op::DirichletBoundaryOperator{MassMatrixMethod}, ctx)

Return a diagonal `Diagonal` mass matrix with zeros at boundary DOFs.
Collected by `build_unsplit_problem` and passed to `ODEFunction`.
"""
function mass_matrix(op::DirichletBoundaryOperator{<:MassMatrixMethod}, ctx::SystemContext)
    layout = ctx.layout
    nx     = ctx.nx
    nvars  = nvariables(layout)
    n      = nvars * nx

    m    = ones(Float64, n)
    vars = op.selector(layout)
    ix   = op.side === :left ? 1 : nx

    @inbounds for ivar in vars
        idx = (ix - 1) * nvars + ivar
        m[idx] = 0.0
    end

    return Diagonal(m)
end


# ---------------------------------------------------------------------------
# CallbackMethod  rhs! / jacobian!  (no-ops; effect delivered via callback)
# ---------------------------------------------------------------------------

function rhs!(du, op::DirichletBoundaryOperator{<:CallbackMethod}, u, ctx::SystemContext, t)
    return du  # no ODE contribution — callback handles enforcement
end

function jacobian!(J, op::DirichletBoundaryOperator{<:CallbackMethod}, u, ctx::SystemContext, t)
    return J   # no Jacobian contribution
end

"""
    build_solver_callback(op::DirichletBoundaryOperator{CallbackMethod}, model) -> DiscreteCallback

Return a `DiscreteCallback` that resets the boundary node to `g(t)` after each
accepted step.  Collected by `solve_problem(model, ...)` and passed to the solver.
"""
function build_solver_callback(op::DirichletBoundaryOperator{<:CallbackMethod},
                                model::SystemModel)
    layout = model.layout
    nx     = model.context.nx
    nvars  = nvariables(layout)
    vars   = op.selector(layout)
    ix     = op.side === :left ? 1 : nx

    idxs = [(ix - 1) * nvars + ivar for ivar in vars]

    condition(u, t, integrator) = true
    function affect!(integrator)
        g = op.bc_fn(integrator.t)
        @inbounds for idx in idxs
            integrator.u[idx] = g
        end
    end

    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end


# ---------------------------------------------------------------------------
# EliminatedMethod  rhs! / jacobian!
# ---------------------------------------------------------------------------

function rhs!(du, op::DirichletBoundaryOperator{<:EliminatedMethod}, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    dx     = ctx.mesh.dx
    invdx2 = inv(dx * dx)

    T_val = op.temperature !== nothing ?
        Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN

    U  = state_view(u, layout, nx)
    dU = state_view(du, layout, nx)

    vars = op.selector(layout)
    g    = op.bc_fn(t)

    # Finite-difference approximation of g'(t)
    eps_t  = 1e-7 * max(1.0, abs(t))
    g_dot  = (op.bc_fn(t + eps_t) - g) / eps_t

    if op.side === :left
        @inbounds for ivar in vars
            D = _eval_D(op.coefficients, ivar, 1, T_val)
            # 1. Cancel the Neumann stencil that LinearDiffusionOperator added at node 1
            dU[ivar, 1] -= D * (U[ivar, 2] - U[ivar, 1]) * invdx2
            # 2. Replace with g'(t) so u[1] tracks the BC exactly
            dU[ivar, 1] += g_dot
            # 3. Correct node 2's stencil to see the exact BC value (not the evolved U[1])
            dU[ivar, 2] += D * (g - U[ivar, 1]) * invdx2
        end
    else  # :right
        @inbounds for ivar in vars
            D = _eval_D(op.coefficients, ivar, 1, T_val)
            dU[ivar, nx] -= D * (U[ivar, nx - 1] - U[ivar, nx]) * invdx2
            dU[ivar, nx] += g_dot
            dU[ivar, nx - 1] += D * (g - U[ivar, nx]) * invdx2
        end
    end

    return du
end

function jacobian!(J, op::DirichletBoundaryOperator{<:EliminatedMethod}, u, ctx::SystemContext, t)
    layout = ctx.layout
    nx     = ctx.nx
    dx     = ctx.mesh.dx
    invdx2 = inv(dx * dx)
    nvars  = nvariables(layout)

    T_val = op.temperature !== nothing ?
        Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN

    vars = op.selector(layout)

    if op.side === :left
        @inbounds for ivar in vars
            D = _eval_D(op.coefficients, ivar, 1, T_val)

            r1 = ivar               # row for node 1
            r2 = nvars + ivar       # row for node 2

            # Cancel the Neumann Jacobian at node 1: LDO added J[1,1]+=-D/dx², J[1,2]+=+D/dx²
            J[r1, r1] += D * invdx2   # cancels the -D/dx² from LDO
            J[r1, r2] -= D * invdx2   # cancels the +D/dx² from LDO

            # Correct node-2 Jacobian: the correction term +D*(g-U[1])/dx² has
            # d/d(U[1]) = -D/dx², which cancels the +D/dx² from the LDO interior stencil.
            J[r2, r1] -= D * invdx2
        end
    else  # :right
        @inbounds for ivar in vars
            D = _eval_D(op.coefficients, ivar, 1, T_val)

            r_nx  = (nx - 1) * nvars + ivar
            r_nm1 = (nx - 2) * nvars + ivar

            J[r_nx, r_nx]  += D * invdx2
            J[r_nx, r_nm1] -= D * invdx2
            J[r_nm1, r_nx] -= D * invdx2
        end
    end

    return J
end


# ===========================================================================
# Generic fallback for build_solver_callback (returns nothing for non-callback operators)
# ===========================================================================

"""
    build_solver_callback(op, model) -> nothing

Default: no callback.  Override for `CallbackMethod` operators.
"""
build_solver_callback(::AbstractOperator, ::SystemModel) = nothing


# ===========================================================================
# Surface flux computation
# ===========================================================================

"""
    surface_fluxes(op, result) -> Dict{Symbol, NamedTuple}

Compute outward diffusive fluxes at both domain boundaries for each diffusing
variable over all saved times.
"""
surface_fluxes(::AbstractDiffusionOperator, result) = Dict{Symbol, NamedTuple}()

function _compute_surface_fluxes(selector, coefficients, temperature,
                                  result::SimulationResult)
    model  = result.model
    layout = model.layout
    mesh   = model.context.mesh
    dx     = mesh.dx
    nx     = model.context.nx
    vars   = selector(layout)
    names  = variable_names(layout)
    nt     = length(result.solution.u)
    ctx    = model.context

    out = Dict{Symbol, NamedTuple}()

    for ivar in vars
        left_flux  = zeros(Float64, nt)
        right_flux = zeros(Float64, nt)

        for it in 1:nt
            t     = result.solution.t[it]
            T_val = temperature !== nothing ?
                Float64(temperature_at(temperature, ctx, t, 1)) : NaN
            D = _eval_D(coefficients, ivar, 1, T_val)

            U = state_view(result.solution.u[it], layout, nx)
            left_flux[it]  = D * (U[ivar, 2]      - U[ivar, 1])  / dx
            right_flux[it] = D * (U[ivar, nx - 1] - U[ivar, nx]) / dx
        end

        out[names[ivar]] = (left = left_flux, right = right_flux)
    end

    return out
end

function surface_fluxes(op::LinearDiffusionOperator, result::SimulationResult)
    return _compute_surface_fluxes(op.selector, op.coefficients, op.temperature, result)
end

function surface_fluxes(op::WeakDirichletBoundaryOperator, result::SimulationResult)
    return _compute_surface_fluxes(op.selector, op.coefficients, op.temperature, result)
end

function surface_fluxes(op::DirichletBoundaryOperator{<:EliminatedMethod},
                         result::SimulationResult)
    return _compute_surface_fluxes(op.selector, op.coefficients, op.temperature, result)
end
