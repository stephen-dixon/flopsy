using OrdinaryDiffEq

function run_simulation(config_path::AbstractString)
    cfg = load_config(config_path)

    # Minimal starter runner using a toy generic reaction-diffusion model.
    length_domain = get(cfg, "length", 1.0)
    nx = get(cfg, "nx", 100)
    tspan = Tuple(get(cfg, "tspan", [0.0, 1.0]))

    mesh = Mesh1D(length_domain, nx)

    vars = [
        VariableInfo(:u, :state, Set([:reaction, :diffusion])),
    ]
    layout = VariableLayout(vars)

    reaction = nothing

    diffusion_coeffs = [get(cfg, "diffusion_coefficient", 1.0)]
    selector(layout) = variables_with_tag(layout, :diffusion)
    diffusion = LinearDiffusionOperator(diffusion_coeffs, selector, nothing)

    model = build_rd_model(
        layout=layout,
        mesh=mesh,
        reaction=reaction,
        diffusion=diffusion,
    )

    n = nvariables(layout) * nx
    u0 = zeros(n)

    # Simple initial pulse
    U0 = state_view(u0, layout, nx)
    mid = cld(nx, 2)
    U0[1, mid] = 1.0

    solver_config = SolverConfig(
        formulation=UnsplitFormulation(),
        algorithm=Rodas5(),
        abstol=1e-8,
        reltol=1e-6,
    )

    return solve_problem(model, u0, tspan, solver_config)
end
