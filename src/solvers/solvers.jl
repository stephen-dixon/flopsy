"""
    solve_problem(model, u0, tspan, solver_config) -> SciML solution
    solve_problem(prob, solver_config) -> SciML solution

Run the time integration.  The first form builds the `ODEProblem` from the model
and then solves it; the second form accepts a pre-built SciML problem.

The first form also collects any `DiscreteCallback`s registered by
`CallbackMethod` Dirichlet operators and merges them into the solve call.
"""
function solve_problem(model::SystemModel, u0, tspan, solver_config::SolverConfig)
    prob = build_problem(model, u0, tspan, solver_config.formulation, solver_config)

    # Collect DiscreteCallbacks from CallbackMethod operators.
    callbacks = _collect_solver_callbacks(model)

    if isempty(callbacks)
        return solve_problem(prob, solver_config)
    else
        cb = CallbackSet(callbacks...)
        kw = merge(solver_config.kwargs, (; callback = cb))
        tmp_config = SolverConfig(
            formulation = solver_config.formulation,
            algorithm = solver_config.algorithm,
            abstol = solver_config.abstol,
            reltol = solver_config.reltol,
            saveat = solver_config.saveat,
            dt = solver_config.dt,
            show_progress = solver_config.show_progress,
            show_solver_stats = solver_config.show_solver_stats,
            write_convergence_trace = solver_config.write_convergence_trace,
            kwargs = kw
        )
        return solve_problem(prob, tmp_config)
    end
end

"""
    _collect_solver_callbacks(model) -> Vector{DiscreteCallback}

Scan all active operators for `CallbackMethod` Dirichlet BCs and collect their
`DiscreteCallback`s.
"""
function _collect_solver_callbacks(model::SystemModel)
    cbs = []
    for op in values(model.operators)
        op isa Nothing && continue
        cb = build_solver_callback(op, model)
        cb !== nothing && push!(cbs, cb)
    end
    return cbs
end
