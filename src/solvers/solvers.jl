function solve_problem(model::SystemModel, u0, tspan, solver_config::SolverConfig)
    prob = build_problem(model, u0, tspan, solver_config.formulation, solver_config)
    return solve_problem(prob, solver_config)
end
