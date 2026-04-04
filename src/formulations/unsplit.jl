function build_problem(model::SystemModel, u0, tspan, ::UnsplitFormulation, solver_config::SolverConfig)
    return build_unsplit_problem(model, u0, tspan, solver_config)
end
