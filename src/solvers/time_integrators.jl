"""
    solve_problem(prob::SciMLBase.AbstractODEProblem, solver_config) -> SciML solution

Run a SciML ODE/DAE problem with the parameters from `solver_config`.
"""
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

"""
    solve_problem(prob::SplitProblem, solver_config) -> SplitSolution

Execute the operator-splitting loop.  Delegates to `_solve_split`.
"""
function solve_problem(prob::SplitProblem, solver_config::SolverConfig)
    formulation = solver_config.formulation
    formulation isa SplitFormulation || throw(ArgumentError(
        "solve_problem(::SplitProblem) requires a SplitFormulation solver config"
    ))

    solver_config.dt !== nothing || throw(ArgumentError(
        "SplitFormulation requires solver_config.dt to be set (macro-step size)"
    ))

    return _solve_split(prob, solver_config, formulation.scheme)
end
