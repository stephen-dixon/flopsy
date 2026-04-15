abstract type AbstractProblemTemplate end

struct Diffusion1DTemplate <: AbstractProblemTemplate end

function instantiate(::Diffusion1DTemplate, cfg::ProblemConfig)
    mesh = _build_mesh(cfg.mesh)

    variables = [VariableInfo(:u, :state, Set([:diffusion, :reaction]))]
    layout = VariableLayout(variables)

    reaction_rate = _get_parameter(cfg, :reaction_rate, 0.0)
    reaction = ToyReactionOperator(reaction_rate)

    coeffs = ConstantDiffusion([_get_parameter(cfg, :diffusion_coefficient, 0.1)])
    selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)
    diffusion = LinearDiffusionOperator(coeffs, selector, nothing)
    boundary = _build_boundary_operator(cfg, layout, coeffs)

    model = build_rd_model(
        layout = layout,
        mesh = mesh,
        reaction = reaction,
        diffusion = diffusion,
        boundary = boundary
    )

    u0 = zeros(Float64, nvariables(layout) * mesh_nx(mesh))
    U0 = state_view(u0, layout, mesh_nx(mesh))
    U0[1, cld(mesh_nx(mesh), 2)] = _get_parameter(cfg, :initial_pulse_amplitude, 1.0)

    return _assemble_problem(cfg, model, u0)
end
