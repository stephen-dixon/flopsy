supports_rhs(::AbstractOperator) = false
supports_implicit_rhs(::AbstractOperator) = false
supports_step(::AbstractOperator) = false
supports_residual(::AbstractOperator) = false
supports_mass_matrix(::AbstractOperator) = false
supports_jacobian(::AbstractOperator) = false

struct NullOperator <: AbstractOperator end

function rhs!(du, ::NullOperator, u, ctx, t)
    return du
end

function implicit_rhs!(du, ::NullOperator, u, ctx, t)
    return du
end

function step!(u, ::NullOperator, ctx, dt, t)
    return u
end

function residual!(res, ::NullOperator, du, u, ctx, t)
    return res
end

mass_matrix(::NullOperator, ctx) = nothing

function jacobian!(J, ::AbstractOperator, u, ctx, t)
    throw(MethodError(jacobian!, (J, typeof(u), typeof(ctx), typeof(t))))
end

function rhs!(du, op::AbstractOperator, u, ctx, t)
    throw(MethodError(rhs!, (du, op, u, ctx, t)))
end

function implicit_rhs!(du, op::AbstractOperator, u, ctx, t)
    throw(MethodError(implicit_rhs!, (du, op, u, ctx, t)))
end

function step!(u, op::AbstractOperator, ctx, dt, t)
    throw(MethodError(step!, (u, op, ctx, dt, t)))
end

function residual!(res, op::AbstractOperator, du, u, ctx, t)
    throw(MethodError(residual!, (res, op, du, u, ctx, t)))
end

mass_matrix(op::AbstractOperator, ctx) = nothing

function active_operators(model::SystemModel)
    return [op for op in values(model.operators) if !(op isa Nothing)]
end
