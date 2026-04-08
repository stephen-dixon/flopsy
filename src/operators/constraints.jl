"""
    ConstraintOperator(f!)

DAE-style constraint operator.  `f!(res, du, u, ctx, t)` writes the algebraic
residual into `res`.  Used with `ResidualFormulation` (not yet fully implemented).
"""
struct ConstraintOperator{F} <: AbstractConstraintOperator
    f!::F
end

supports_residual(::ConstraintOperator) = true

function residual!(res, op::ConstraintOperator, du, u, ctx, t)
    op.f!(res, du, u, ctx, t)
    return res
end
