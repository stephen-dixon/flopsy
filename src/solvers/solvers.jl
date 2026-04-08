"""
    solve_problem(model, u0, tspan, solver_config) -> SciML solution
    solve_problem(prob, solver_config) -> SciML solution

Run the time integration.  The first form builds the `ODEProblem` from the model
and then solves it; the second form accepts a pre-built SciML problem.
"""
function solve_problem(model::SystemModel, u0, tspan, solver_config::SolverConfig)
    prob = build_problem(model, u0, tspan, solver_config.formulation, solver_config)
    return solve_problem(prob, solver_config)
end
