supports_rhs(::AbstractOperator) = false
supports_implicit_rhs(::AbstractOperator) = false
supports_step(::AbstractOperator) = false
supports_residual(::AbstractOperator) = false
supports_mass_matrix(::AbstractOperator) = false
supports_jacobian(::AbstractOperator) = false

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

function rhs!(du, op::AbstractOperator, u, ctx, t)
    throw(ArgumentError("rhs! not implemented for operator type $(typeof(op))"))
end

function implicit_rhs!(du, op::AbstractOperator, u, ctx, t)
    throw(ArgumentError("implicit_rhs! not implemented for operator type $(typeof(op))"))
end

function step!(u, op::AbstractOperator, ctx, dt, t)
    throw(ArgumentError("step! not implemented for operator type $(typeof(op))"))
end

function residual!(res, op::AbstractOperator, du, u, ctx, t)
    throw(ArgumentError("residual! not implemented for operator type $(typeof(op))"))
end

function jacobian!(J, op::AbstractOperator, u, ctx, t)
    throw(ArgumentError("jacobian! not implemented for operator type $(typeof(op))"))
end

mass_matrix(op::AbstractOperator, ctx) = nothing

function active_operators(model::SystemModel)
    return [op for op in values(model.operators) if !(op isa Nothing)]
end
