"Return `true` if the operator contributes an explicit RHS term `du += f(u,t)`."
supports_rhs(::AbstractOperator) = false

"Return `true` if the operator contributes an implicit RHS term."
supports_implicit_rhs(::AbstractOperator) = false

"Return `true` if the operator implements a custom `step!` update."
supports_step(::AbstractOperator) = false

"Return `true` if the operator can form a DAE residual."
supports_residual(::AbstractOperator) = false

"Return `true` if the operator provides a mass matrix."
supports_mass_matrix(::AbstractOperator) = false

"Return `true` if the operator can compute an analytic Jacobian."
supports_jacobian(::AbstractOperator) = false

"""
    NullOperator

A no-op operator that supports all capability flags.  Used as a placeholder when
an operator slot (reaction, diffusion, etc.) is not needed.
"""
struct NullOperator <: AbstractOperator end

supports_rhs(::NullOperator) = true
supports_implicit_rhs(::NullOperator) = true
supports_step(::NullOperator) = true
supports_residual(::NullOperator) = true
supports_mass_matrix(::NullOperator) = true
supports_jacobian(::NullOperator) = true

rhs!(du, ::NullOperator, u, ctx, t) = du
implicit_rhs!(du, ::NullOperator, u, ctx, t) = du
step!(u, ::NullOperator, ctx, dt, t) = u
residual!(res, ::NullOperator, du, u, ctx, t) = res
jacobian!(J, ::NullOperator, u, ctx, t) = J
mass_matrix(::NullOperator, ctx) = nothing
jacobian_node_sparsity(::NullOperator, layout) = Set{Tuple{Int,Int}}()

"""
    rhs!(du, op, u, ctx, t) -> du

Accumulate the explicit right-hand-side contribution of operator `op` into `du`.
`du` is zeroed by `OperatorSum` before each call; individual operators should
add their contribution rather than overwrite.
"""
function rhs!(du, op::AbstractOperator, u, ctx, t)
    throw(ArgumentError("rhs! not implemented for operator type $(typeof(op))"))
end

"""
    implicit_rhs!(du, op, u, ctx, t) -> du

Accumulate the implicit RHS contribution of operator `op` into `du`.
"""
function implicit_rhs!(du, op::AbstractOperator, u, ctx, t)
    throw(ArgumentError("implicit_rhs! not implemented for operator type $(typeof(op))"))
end

"""
    step!(u, op, ctx, dt, t) -> u

Apply a custom operator-specific update step of size `dt` in-place to `u`.
"""
function step!(u, op::AbstractOperator, ctx, dt, t)
    throw(ArgumentError("step! not implemented for operator type $(typeof(op))"))
end

"""
    residual!(res, op, du, u, ctx, t) -> res

Evaluate the DAE residual `F(du, u, t)` for operator `op` into `res`.
"""
function residual!(res, op::AbstractOperator, du, u, ctx, t)
    throw(ArgumentError("residual! not implemented for operator type $(typeof(op))"))
end

"""
    jacobian!(J, op, u, ctx, t) -> J

Fill the Jacobian matrix `J` with ∂f/∂u for operator `op`.
"""
function jacobian!(J, op::AbstractOperator, u, ctx, t)
    throw(ArgumentError("jacobian! not implemented for operator type $(typeof(op))"))
end

"""
    mass_matrix(op, ctx)

Return the mass matrix for operator `op`, or `nothing` if not applicable.
"""
mass_matrix(op::AbstractOperator, ctx) = nothing

"""
    jacobian_node_sparsity(op, layout) -> Set{Tuple{Int,Int}} or nothing

Return the per-node Jacobian nonzero pattern as a set of `(row, col)` pairs within
the `nvars × nvars` block for a single spatial node.  Return `nothing` to signal a
fully dense block (the safe fallback).

Called by `_build_jac_prototype` to assemble the sparse Jacobian prototype.
Off-node (inter-node diffusion) entries are handled separately via
`diffusion_variable_indices` and do not need to be included here.
Pairs use 1-based variable indices `1..nvariables(layout)`.
"""
jacobian_node_sparsity(::AbstractOperator, layout) = nothing

"""
    active_operators(model) -> Vector

Return all non-`nothing` operator values from `model.operators`.
"""
function active_operators(model::SystemModel)
    return [op for op in values(model.operators) if !(op isa Nothing)]
end
