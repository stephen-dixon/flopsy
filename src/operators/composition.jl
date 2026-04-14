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

# Mass matrix: true if any sub-op has a mass matrix; combined matrix takes elementwise min.
function supports_mass_matrix(op::OperatorSum)
    return any(supports_mass_matrix, op.ops)
end

function mass_matrix(op::OperatorSum, ctx::SystemContext)
    mm_ops = [sub for sub in op.ops if supports_mass_matrix(sub)]
    isempty(mm_ops) && return nothing

    nvars = nvariables(ctx.layout)
    n     = nvars * ctx.nx
    m     = ones(Float64, n)

    for sub in mm_ops
        M_sub = mass_matrix(sub, ctx)
        M_sub === nothing && continue
        @inbounds for i in 1:n
            m[i] = min(m[i], M_sub.diag[i])
        end
    end

    return Diagonal(m)
end

# Callback collection: recurse into sub-operators.
function build_solver_callback(op::OperatorSum, model::SystemModel)
    cbs = []
    for sub in op.ops
        cb = build_solver_callback(sub, model)
        cb !== nothing && push!(cbs, cb)
    end
    isempty(cbs) && return nothing
    length(cbs) == 1 && return cbs[1]
    return CallbackSet(cbs...)
end
