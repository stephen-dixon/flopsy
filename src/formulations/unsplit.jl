function build_problem(model::SystemModel, u0, tspan, formulation::UnsplitFormulation, solver_config)
    return build_unsplit_problem(model, u0, tspan, solver_config)
end
