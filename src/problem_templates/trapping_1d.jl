struct Trapping1DTemplate <: AbstractProblemTemplate end

function instantiate(::Trapping1DTemplate, cfg::ProblemConfig)
    mesh = _build_mesh(cfg.mesh)
    diffusion = ConstantDiffusion([
        _get_parameter(cfg, :diffusion_coefficient, 0.01),
        0.0,
    ])
    model = build_trapping_model(
        mesh = mesh,
        k_trap = _get_parameter(cfg, :k_trap, 5.0),
        k_detrap = _get_parameter(cfg, :k_detrap, 0.5),
        diffusion_coefficient = diffusion,
    )

    u0 = zeros(Float64, nvariables(model.layout) * mesh_nx(mesh))
    U0 = state_view(u0, model.layout, mesh_nx(mesh))
    U0[1, cld(mesh_nx(mesh), 2)] = _get_parameter(cfg, :initial_mobile_pulse_amplitude, 1.0)
    U0[2, :] .= _get_parameter(cfg, :initial_trap_occupancy, 0.0)

    return _assemble_problem(cfg, model, u0)
end
