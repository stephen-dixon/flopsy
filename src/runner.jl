module Runner

using ..Config
using ..Mesh
using ..Types
using ..State
using ..Solver
using ..Output

export run_simulation

function run_simulation(config_path::String)
    config = Config.load_config(config_path)

    x, dx = Mesh.create_mesh(config.nx, config.length)

    species = Types.SpeciesInfo(["mobile"], 1)
    ctx = Types.SimulationContext(species, config.nx, dx)

    u0 = State.create_initial_state(length(species.names), species.nlevels, config.nx)

    sol = Solver.solve_problem(u0, config.tspan, ctx, config)

    Output.write_outputs(sol, ctx, config)

    return sol
end

end
