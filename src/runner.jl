"""
    run_simulation(config_path) -> SimulationResult

Load a TOML config file and run the simulation described by it.
Returns a `SimulationResult` wrapping the solution.
"""
function run_simulation(config_path::AbstractString)
    cfg = load_config(config_path)

    model_type = get(cfg, "model_type", "toy_rd")

    length_domain = get(cfg, "length", 1.0)
    nx = get(cfg, "nx", 101)
    tspan_vec = get(cfg, "tspan", [0.0, 1.0])
    tspan = (tspan_vec[1], tspan_vec[2])

    mesh = Mesh1D(length_domain, nx)

    model = if model_type == "toy_rd"
        build_toy_rd_model(cfg, mesh)
    elseif model_type == "toy_trapping"
        build_toy_trapping_model(cfg, mesh)
    elseif model_type == "fake_hotgates_trapping"
        build_fake_hotgates_trapping_model(cfg, mesh)
    else
        throw(ArgumentError("Unknown model_type: $model_type"))
    end

    u0 = build_initial_state(cfg, model)
    algorithm = _build_algorithm(cfg)

    solver_config = SolverConfig(
        formulation = UnsplitFormulation(),
        algorithm = algorithm,
        abstol = get(cfg, "abstol", 1e-8),
        reltol = get(cfg, "reltol", 1e-6),
        saveat = get(cfg, "saveat", nothing),
    )

    print_run_banner(cfg, solver_config, model)

    t_start = time()
    sol = solve_problem(model, u0, tspan, solver_config)
    elapsed_time = time() - t_start

    metadata = Dict{String,Any}(
        "model_type" => model_type,
        "algorithm" => string(typeof(algorithm)),
        "formulation" => string(typeof(solver_config.formulation)),
        "nx" => model.context.nx,
        "nvariables" => nvariables(model.layout),
        "retcode" => string(sol.retcode),
        "elapsed_walltime_s" => elapsed_time
    )

    summaries = Dict{Symbol,Any}()

    # mock example of an extra derived timeseries for trapping-like runs
    if nvariables(model.layout) >= 2
        summaries[:extra_timeseries] = Dict(
            "integral_total_inventory" =>
                integrated_variable(wrap_result(model, sol, cfg), 1) .+
                integrated_variable(wrap_result(model, sol, cfg), 2)
        )
    end

    return wrap_result(model, sol, cfg; summaries=summaries, metadata=metadata)
end


function build_toy_rd_model(cfg, mesh::Mesh1D)
    diffusion_coefficient = get(cfg, "diffusion_coefficient", 0.1)
    reaction_rate = get(cfg, "reaction_rate", 1.0)

    variables = [
        VariableInfo(:u, :state, Set([:reaction, :diffusion])),
    ]
    layout = VariableLayout(variables)

    reaction = ToyReactionOperator(reaction_rate)

    selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)
    diffusion = LinearDiffusionOperator([diffusion_coefficient], selector, nothing)

    return build_rd_model(
        layout = layout,
        mesh = mesh,
        reaction = reaction,
        diffusion = diffusion,
    )
end


function build_toy_trapping_model(cfg, mesh::Mesh1D)
    k_trap = get(cfg, "k_trap", 5.0)
    k_detrap = get(cfg, "k_detrap", 0.5)
    diffusion_coefficient = get(cfg, "diffusion_coefficient", 0.01)

    return build_trapping_model(
        mesh = mesh,
        k_trap = k_trap,
        k_detrap = k_detrap,
        diffusion_coefficient = diffusion_coefficient,
        mobile_name = :c,
        trap_name = :theta,
    )
end


function build_initial_state(cfg, model::SystemModel)
    layout = model.layout
    nx = model.context.nx

    u0 = zeros(nvariables(layout) * nx)
    U0 = state_view(u0, layout, nx)

    if nvariables(layout) == 1
        mid = cld(nx, 2)
        U0[1, mid] = get(cfg, "initial_pulse_amplitude", 1.0)
    elseif nvariables(layout) == 2
        mid = cld(nx, 2)
        U0[1, mid] = get(cfg, "initial_mobile_pulse_amplitude", 1.0)
        U0[2, :] .= get(cfg, "initial_trap_occupancy", 0.0)
    else
        throw(ArgumentError("No default initial-state builder for $(nvariables(layout)) variables"))
    end

    return u0
end


function _build_algorithm(cfg)
    alg = get(cfg, "algorithm", "Rodas5")
    # alg = get(cfg, "algorithm", "Rodas5")

    if alg == "Rodas5"
        return Rodas5(autodiff = AutoFiniteDiff())
    elseif alg == "CVODE_BDF"
        return Sundials.CVODE_BDF()
    else
        throw(ArgumentError("Unknown algorithm: $alg"))
    end
end


function build_fake_hotgates_trapping_model(cfg, mesh::Mesh1D)
    k_trap = get(cfg, "k_trap", 5.0)
    k_detrap = get(cfg, "k_detrap", 0.5)
    diffusion_coefficient = get(cfg, "diffusion_coefficient", 0.01)
    temperature = get(cfg, "temperature", 300.0)

    backend = FakeHotgatesModel(k_trap, k_detrap)

    nx = length(mesh.x)

    adaptor = HotgatesTrappingAdaptor(
        [1],               # mobile indices in local Flopsy state
        [2],               # trap indices in local Flopsy state
        ["c"],
        ["theta"],
        String[],          # no defect species
        zeros(Float64, 0, nx),  # FakeHotgatesModel does not use defects
    )

    diffusion_coeffs = [diffusion_coefficient, 0.0]

    return build_hotgates_trapping_model(
        mesh = mesh,
        model = backend,
        adaptor = adaptor,
        temperature = ConstantTemperature(temperature),
        diffusion_coefficients = diffusion_coeffs,
    )
end
