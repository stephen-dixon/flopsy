"""
    build_problem(model, u0, tspan, formulation, solver_config) -> SciML problem

Assemble a SciML problem (e.g. `ODEProblem`) from the model, initial state,
time span, and formulation.  Dispatches on the formulation type.
"""
function build_problem(model::SystemModel, u0, tspan, ::UnsplitFormulation, solver_config::SolverConfig)
    return build_unsplit_problem(model, u0, tspan, solver_config)
end

# IMEXFormulation and IMEXReactionFormulation are dispatched in formulations/imex.jl

function build_unsplit_problem(model::SystemModel, u0, tspan, solver_config::SolverConfig)
    ops = Tuple(active_operators(model))
    total_op = OperatorSum(ops)

    supports_rhs(total_op) || throw(ArgumentError(
        "Unsplit formulation requires rhs! support for all active operators"
    ))

    function f!(du, u, p, t)
        rhs!(du, total_op, u, model.context, t)
        return nothing
    end

    # Collect mass matrix from MassMatrixMethod operators (if any).
    M = _collect_mass_matrix(ops, model.context)

    if supports_jacobian(total_op)
        prototype = _build_jac_prototype(model, ops)

        function jac!(J, u, p, t)
            jacobian!(J, total_op, u, model.context, t)
            return nothing
        end

        if M !== nothing
            ode_f = ODEFunction(f!, jac = jac!, jac_prototype = prototype,
                                mass_matrix = M)
        else
            ode_f = ODEFunction(f!, jac = jac!, jac_prototype = prototype)
        end
    else
        if M !== nothing
            ode_f = ODEFunction(f!, mass_matrix = M)
        else
            ode_f = ODEFunction(f!)
        end
    end

    return ODEProblem(ode_f, u0, tspan)
end


# ---------------------------------------------------------------------------
# Mass matrix collection
# ---------------------------------------------------------------------------

"""
    _collect_mass_matrix(ops, ctx) -> Diagonal or nothing

Scan operators for `MassMatrixMethod` Dirichlet BCs.  Builds a combined
diagonal mass matrix (identity everywhere except boundary DOFs which get 0).
Returns `nothing` if no mass-matrix operators are found.
"""
function _collect_mass_matrix(ops, ctx::SystemContext)
    any(supports_mass_matrix, ops) || return nothing

    layout = ctx.layout
    nx     = ctx.nx
    nvars  = nvariables(layout)
    n      = nvars * nx

    m = ones(Float64, n)
    for op in ops
        supports_mass_matrix(op) || continue
        M_op = mass_matrix(op, ctx)
        M_op === nothing && continue
        @inbounds for i in 1:n
            m[i] = min(m[i], M_op.diag[i])
        end
    end

    return Diagonal(m)
end


# ---------------------------------------------------------------------------
# Sparse Jacobian prototype
# ---------------------------------------------------------------------------

function _build_jac_prototype(model::SystemModel, ops)
    layout = model.context.layout
    nx     = model.context.nx
    nvars  = nvariables(layout)
    n      = nvars * nx

    entries = Set{Tuple{Int,Int}}()

    for ix in 1:nx
        offset = (ix - 1) * nvars
        for iv1 in 1:nvars, iv2 in 1:nvars
            push!(entries, (offset + iv1, offset + iv2))
        end
    end

    for op in ops
        for ivar in diffusion_variable_indices(op, layout)
            for ix in 1:(nx - 1)
                r = (ix - 1) * nvars + ivar
                c = ix       * nvars + ivar
                push!(entries, (r, c))
                push!(entries, (c, r))
            end
        end
    end

    Is = [e[1] for e in entries]
    Js = [e[2] for e in entries]
    return sparse(Is, Js, ones(Float64, length(Is)), n, n)
end
