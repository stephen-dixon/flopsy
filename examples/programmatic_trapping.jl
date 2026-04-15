using Flopsy
using OrdinaryDiffEq

mesh = Mesh1D(1.0, 51)
model = build_trapping_model(
    mesh = mesh,
    k_trap = 2.0,
    k_detrap = 0.25,
    diffusion_coefficient = 0.05
)

u0 = zeros(nvariables(model.layout) * length(mesh.x))
tspan = (0.0, 0.25)

config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm = Rodas5P(),
    abstol = 1e-8,
    reltol = 1e-6,
    saveat = [0.0, 0.125, 0.25]
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
println("retcode: ", getproperty(result.solution, :retcode))
