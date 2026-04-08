"""
    build_problem(model, u0, tspan, formulation, solver_config) -> SciML problem

Assemble a SciML problem (e.g. `ODEProblem`) from the model, initial state,
time span, and formulation.  Dispatches on the formulation type.
"""
function build_problem(model::SystemModel, u0, tspan, ::UnsplitFormulation, solver_config::SolverConfig)
    return build_unsplit_problem(model, u0, tspan, solver_config)
end

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

    if supports_jacobian(total_op)
        prototype = _build_jac_prototype(model, ops)

        function jac!(J, u, p, t)
            jacobian!(J, total_op, u, model.context, t)
            return nothing
        end

        ode_f = ODEFunction(f!, jac = jac!, jac_prototype = prototype)
    else
        ode_f = ODEFunction(f!)
    end

    return ODEProblem(ode_f, u0, tspan)
end


# ---------------------------------------------------------------------------
# Sparse Jacobian prototype
# ---------------------------------------------------------------------------

# Build a sparse Float64 matrix whose sparsity pattern reflects the known
# structure of the assembled operator:
# - Within each node: all-to-all (dense nvars×nvars block) for reaction coupling.
# - Between adjacent nodes: same-variable entries for diffusion coupling.
function _build_jac_prototype(model::SystemModel, ops)
    layout = model.context.layout
    nx     = model.context.nx
    nvars  = nvariables(layout)
    n      = nvars * nx

    entries = Set{Tuple{Int,Int}}()

    # Dense within-node blocks: every variable at node ix can couple to every
    # other variable at node ix (reaction terms).
    for ix in 1:nx
        offset = (ix - 1) * nvars
        for iv1 in 1:nvars, iv2 in 1:nvars
            push!(entries, (offset + iv1, offset + iv2))
        end
    end

    # Spatial coupling: same variable at adjacent nodes (diffusion operators).
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
