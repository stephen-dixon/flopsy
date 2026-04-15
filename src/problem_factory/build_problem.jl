"""
    build_problem(cfg::ProblemConfig) -> SimulationProblem

Build a fully assembled `SimulationProblem` from a typed configuration.
"""
function build_problem(cfg::ProblemConfig)
    validate(cfg)
    template = _select_problem_template(cfg.problem_type)
    return instantiate(template, cfg)
end

function _select_problem_template(problem_type::Symbol)
    if problem_type == :diffusion_1d
        return Diffusion1DTemplate()
    elseif problem_type == :trapping_1d
        return Trapping1DTemplate()
    elseif problem_type == :hotgates_trapping
        return HotgatesTrappingTemplate()
    end

    throw(ArgumentError("No problem template registered for $(problem_type)"))
end

function _assemble_problem(cfg::ProblemConfig, model::SystemModel, u0::Vector{Float64})
    tspan = (_get_parameter(cfg, :t0, 0.0), _get_parameter(cfg, :tend, 1.0))
    solver_config = _build_solver_config(cfg.solver)
    return SimulationProblem(model, u0, tspan, solver_config, cfg)
end

function _build_solver_config(cfg::InputSolverConfig)
    return SolverConfig(
        formulation = _build_formulation(cfg.formulation),
        algorithm = _build_algorithm(cfg.algorithm),
        abstol = cfg.abstol,
        reltol = cfg.reltol,
        saveat = cfg.saveat,
        dt = cfg.dt,
        show_progress = false,
        show_solver_stats = false
    )
end

function _build_formulation(formulation::Symbol)
    if formulation == :unsplit
        return UnsplitFormulation()
    elseif formulation == :imex
        return IMEXFormulation()
    elseif formulation == :imex_reaction
        return IMEXReactionFormulation()
    elseif formulation == :split
        return SplitFormulation(StrangSplit())
    elseif formulation == :residual
        return ResidualFormulation()
    end

    throw(ArgumentError("Unsupported solver formulation $(formulation)"))
end

function _build_algorithm(algorithm::Symbol)
    if algorithm == :Rodas5
        return Rodas5(autodiff = AutoFiniteDiff())
    elseif algorithm == :Rodas5P
        return Rodas5P()
    elseif algorithm == :KenCarp4
        return KenCarp4()
    elseif algorithm == :CVODE_BDF
        return Sundials.CVODE_BDF()
    end

    throw(ArgumentError("Unsupported solver algorithm $(algorithm)"))
end

"""
    solve(problem::SimulationProblem) -> SimulationResult

Solve a config-built `SimulationProblem` and return a wrapped result.
"""
function solve(problem::SimulationProblem)
    sol = solve_problem(problem.model, problem.u0, problem.tspan, problem.solver_config)
    metadata = Dict{String, Any}(
        "problem_type" => string(problem.config.problem_type),
        "algorithm" => string(typeof(problem.solver_config.algorithm)),
        "formulation" => string(typeof(problem.solver_config.formulation)),
        "nx" => problem.model.context.nx,
        "nvariables" => nvariables(problem.model.layout),
        "retcode" => string(getproperty(sol, :retcode))
    )
    return wrap_result(problem.model, sol, problem.config; metadata = metadata)
end

"""
    run_simulation(config_path) -> SimulationResult

Load a TOML input deck, build a `SimulationProblem`, and solve it.
"""
function run_simulation(config_path::AbstractString)
    cfg = load_config(config_path)
    problem = build_problem(cfg)
    return solve(problem)
end

function _build_mesh(cfg::MeshConfig)
    return Mesh1D(cfg.xmax - cfg.xmin, cfg.nx)
end

mesh_nx(mesh::Mesh1D) = length(mesh.x)

function _get_parameter(cfg::ProblemConfig, key::Symbol, default::Float64)
    return hasproperty(cfg.parameters, key) ? getproperty(cfg.parameters, key) : default
end

function _build_boundary_operator(cfg::ProblemConfig, layout::VariableLayout, coefficients)
    isempty(cfg.boundary_conditions) && return nothing

    selector = let bcs = cfg.boundary_conditions
        function (layout_inner::VariableLayout)
            vars = Set{Int}()
            for bc in bcs
                idx = findfirst(==(bc.variable), variable_names(layout_inner))
                idx === nothing && continue
                push!(vars, idx)
            end
            return sort!(collect(vars))
        end
    end

    ops = AbstractOperator[]
    for bc in cfg.boundary_conditions
        bc.kind == :dirichlet || continue
        bc_fn = let value = bc.value
            t -> value
        end
        if bc.method == :weak
            if bc.side == :left
                push!(ops, WeakDirichletBoundaryOperator(selector, coefficients, nothing; left = bc_fn))
            else
                push!(ops, WeakDirichletBoundaryOperator(selector, coefficients, nothing; right = bc_fn))
            end
        else
            method = _build_boundary_method(bc.method)
            push!(ops,
                DirichletBoundaryOperator(
                    bc.side, bc_fn, selector, coefficients, nothing; method = method))
        end
    end

    isempty(ops) && return nothing
    length(ops) == 1 && return ops[1]
    return OperatorSum(Tuple(ops))
end

function _build_boundary_method(method::Symbol)
    if method == :penalty
        return PenaltyMethod()
    elseif method == :mass_matrix
        return MassMatrixMethod()
    elseif method == :callback
        return CallbackMethod()
    elseif method == :eliminated
        return EliminatedMethod()
    end
    throw(ArgumentError("Unsupported boundary condition method $(method)"))
end

function _apply_diffusion_ic!(U0, cfg::ProblemConfig)
    ic = cfg.initial_conditions
    nx = size(U0, 2)

    if ic.kind == :uniform
        U0[1, :] .= something(ic.value, _get_parameter(cfg, :initial_uniform_value, 0.0))
    else
        U0[1, cld(nx, 2)] = something(ic.amplitude, _get_parameter(cfg, :initial_pulse_amplitude, 1.0))
    end

    return U0
end

function _apply_trapping_ic!(U0, cfg::ProblemConfig)
    ic = cfg.initial_conditions
    nx = size(U0, 2)
    U0[2, :] .= something(ic.trap_occupancy, _get_parameter(cfg, :initial_trap_occupancy, 0.0))

    if ic.kind == :uniform
        U0[1, :] .= something(ic.mobile_value, _get_parameter(cfg, :initial_mobile_uniform_value, 0.0))
    else
        U0[1, cld(nx, 2)] = something(ic.mobile_amplitude, _get_parameter(cfg, :initial_mobile_pulse_amplitude, 1.0))
    end

    return U0
end
