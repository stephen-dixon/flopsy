"""
    OperatorSum(ops)

Composite operator that accumulates contributions from a tuple of sub-operators.

For each capability (`rhs!`, `jacobian!`, etc.) only sub-operators for which the
corresponding `supports_*` predicate returns `true` are called.  `supports_jacobian`
returns `true` only when **all** sub-operators support it, ensuring the assembled
Jacobian is complete.

`OperatorSum` is used internally by `build_unsplit_problem` to combine a model's
active operators into a single callable.  Users rarely need to construct it directly.
"""
struct OperatorSum{Ops} <: AbstractOperator
    ops::Ops
end

supports_rhs(op::OperatorSum) = all(supports_rhs, op.ops)
supports_implicit_rhs(op::OperatorSum) = all(supports_implicit_rhs, op.ops)
supports_step(op::OperatorSum) = all(supports_step, op.ops)
supports_residual(op::OperatorSum) = all(supports_residual, op.ops)
supports_jacobian(op::OperatorSum) = all(supports_jacobian, op.ops)

function rhs!(du, op::OperatorSum, u, ctx, t)
    fill!(du, zero(eltype(du)))

    tmp = get!(ctx.scratch, :rhs_tmp) do
        similar(du)
    end

    for subop in op.ops
        supports_rhs(subop) || continue
        fill!(tmp, zero(eltype(tmp)))
        rhs!(tmp, subop, u, ctx, t)
        @. du += tmp
    end

    return du
end

function implicit_rhs!(du, op::OperatorSum, u, ctx, t)
    fill!(du, zero(eltype(du)))

    tmp = get!(ctx.scratch, :implicit_rhs_tmp) do
        similar(du)
    end

    for subop in op.ops
        supports_implicit_rhs(subop) || continue
        fill!(tmp, zero(eltype(tmp)))
        implicit_rhs!(tmp, subop, u, ctx, t)
        @. du += tmp
    end

    return du
end

function jacobian!(J, op::OperatorSum, u, ctx, t)
    fill!(J, zero(eltype(J)))
    for subop in op.ops
        jacobian!(J, subop, u, ctx, t)
    end
    return J
end
