struct ConstraintOperator{F} <: AbstractConstraintOperator
    f!::F
end

supports_residual(::ConstraintOperator) = true

function residual!(res, op::ConstraintOperator, du, u, ctx, t)
    op.f!(res, du, u, ctx, t)
    return res
end
