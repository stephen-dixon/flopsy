module Solver

using DifferentialEquations
using Sundials
using ..RHS

export solve_problem

function choose_solver(name::String)
    if name == "Rodas5"
        return Rodas5()
    elseif name == "CVODE_BDF"
        return CVODE_BDF()
    else
        error("Unknown solver: $name")
    end
end

function solve_problem(u0, tspan, ctx, config)
    prob = ODEProblem(RHS.rhs!, u0, tspan, ctx)

    solver = choose_solver(config.solver)

    return solve(prob, solver;
        abstol=config.abstol,
        reltol=config.reltol
    )
end

end
