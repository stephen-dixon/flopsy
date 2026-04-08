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

    return ODEProblem(f!, u0, tspan)
end
