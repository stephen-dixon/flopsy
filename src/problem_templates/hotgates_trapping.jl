struct HotgatesTrappingTemplate <: AbstractProblemTemplate end

function instantiate(::HotgatesTrappingTemplate, cfg::ProblemConfig)
    mesh = _build_mesh(cfg.mesh)
    backend = FakeHotgatesModel(
        _get_parameter(cfg, :k_trap, 5.0),
        _get_parameter(cfg, :k_detrap, 0.5),
    )

    nx = mesh_nx(mesh)
    adaptor = HotgatesTrappingAdaptor(
        [1],
        [2],
        ["c"],
        ["theta"],
        String[],
        zeros(Float64, 0, nx),
    )

    coeffs = [_get_parameter(cfg, :diffusion_coefficient, 0.01), 0.0]
    model = build_hotgates_trapping_model(
        mesh = mesh,
        model = backend,
        adaptor = adaptor,
        temperature = ConstantTemperature(_get_parameter(cfg, :temperature, 300.0)),
        diffusion_coefficients = coeffs,
    )

    u0 = zeros(Float64, nvariables(model.layout) * nx)
    U0 = state_view(u0, model.layout, nx)
    U0[1, cld(nx, 2)] = _get_parameter(cfg, :initial_mobile_pulse_amplitude, 1.0)
    U0[2, :] .= _get_parameter(cfg, :initial_trap_occupancy, 0.0)

    return _assemble_problem(cfg, model, u0)
end
