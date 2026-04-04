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
