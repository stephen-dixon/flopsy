function solve_problem(prob::SciMLBase.AbstractODEProblem, solver_config::SolverConfig)
    kwargs = copy(solver_config.kwargs)

    solver_config.saveat !== nothing && (kwargs[:saveat] = solver_config.saveat)
    solver_config.dt !== nothing && (kwargs[:dt] = solver_config.dt)

    return SciMLBase.solve(
        prob,
        solver_config.algorithm;
        abstol = solver_config.abstol,
        reltol = solver_config.reltol,
        kwargs...,
    )
end
