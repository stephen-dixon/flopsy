using TOML

"""
    load_config(path) -> ProblemConfig

Deprecated legacy entry point. Load a legacy TOML config into a typed
`ProblemConfig`.
"""
function load_config(path::AbstractString)
    Base.depwarn(
        "`load_config` is deprecated. The supported TOML workflow is `parse_input_deck`/`validate_input_deck`/`run_input_deck` with the registry-driven input-deck system.",
        :load_config
    )
    raw = TOML.parsefile(path)
    return parse_config(raw)
end

"""
    parse_config(raw) -> ProblemConfig

Deprecated legacy entry point. Parse an `AbstractDict` into the typed config
objects used by the old problem factory.
"""
function parse_config(raw::AbstractDict)
    Base.depwarn(
        "`parse_config` is deprecated. The supported TOML workflow is the registry-driven input-deck system via `parse_input_deck`.",
        :parse_config
    )
    cfg = _normalize_legacy_config(raw)
    mesh = _parse_mesh_config(_required_table(cfg, "mesh"))
    solver = _parse_input_solver_config(_required_table(cfg, "solver"))
    bcs = _parse_boundary_conditions(get(cfg, "boundary_conditions", Any[]))
    initial_conditions = _parse_initial_conditions(get(cfg, "initial_conditions", Dict{
        String, Any}()))
    parameters = _parse_parameters(cfg)

    return ProblemConfig(
        Symbol(_required_string(cfg, "problem_type")),
        mesh,
        solver,
        bcs,
        initial_conditions,
        parameters
    )
end

function _normalize_legacy_config(raw::AbstractDict)
    haskey(raw, "problem_type") && return raw
    haskey(raw, "model_type") || return raw

    problem_type = get(raw, "model_type", "diffusion_1d")
    tspan_raw = get(raw, "tspan", [0.0, 1.0])
    length_domain = get(raw, "length", 1.0)
    saveat_raw = get(raw, "saveat", nothing)

    formulation = get(raw, "formulation", "unsplit")
    algorithm = get(raw, "algorithm", "Rodas5")

    normalized = Dict{String, Any}()
    normalized["problem_type"] = _normalize_problem_type(problem_type)
    normalized["mesh"] = Dict(
        "kind" => "uniform_1d",
        "nx" => get(raw, "nx", 101),
        "xmin" => 0.0,
        "xmax" => length_domain
    )
    normalized["solver"] = Dict(
        "formulation" => formulation,
        "algorithm" => algorithm,
        "dt" => get(raw, "dt", nothing),
        "abstol" => get(raw, "abstol", 1e-8),
        "reltol" => get(raw, "reltol", 1e-6),
        "saveat" => saveat_raw
    )
    normalized["parameters"] = Dict(
        "t0" => tspan_raw[1],
        "tend" => tspan_raw[2]
    )
    normalized["initial_conditions"] = Dict{String, Any}()

    passthrough = (
        "diffusion_coefficient",
        "reaction_rate",
        "k_trap",
        "k_detrap",
        "temperature",
        "initial_pulse_amplitude",
        "initial_mobile_pulse_amplitude",
        "initial_trap_occupancy"
    )
    for key in passthrough
        if haskey(raw, key)
            normalized["parameters"][key] = raw[key]
        end
    end

    return normalized
end

function _parse_mesh_config(raw::AbstractDict)
    return MeshConfig(
        Symbol(get(raw, "kind", "uniform_1d")),
        Int(get(raw, "nx", 101)),
        Float64(get(raw, "xmin", 0.0)),
        Float64(get(raw, "xmax", get(raw, "length", 1.0)))
    )
end

function _parse_input_solver_config(raw::AbstractDict)
    saveat = _parse_saveat(get(raw, "saveat", nothing))
    dt_value = get(raw, "dt", nothing)
    dt = dt_value === nothing ? nothing : Float64(dt_value)

    return InputSolverConfig(
        _parse_formulation_symbol(get(raw, "formulation", get(raw, "method", "unsplit"))),
        _parse_algorithm_spec(get(raw, "algorithm", get(raw, "method", "Rodas5"))),
        dt,
        Float64(get(raw, "abstol", 1e-8)),
        Float64(get(raw, "reltol", 1e-6)),
        saveat
    )
end

function _parse_boundary_conditions(raw)
    isa(raw, AbstractVector) ||
        throw(ArgumentError("boundary_conditions must be an array of tables"))
    bcs = BoundaryConditionConfig{Float64}[]
    for item in raw
        item isa AbstractDict ||
            throw(ArgumentError("boundary condition entries must be TOML tables"))
        push!(bcs,
            BoundaryConditionConfig(
                Symbol(get(item, "variable", "u")),
                Symbol(_required_string(item, "side")),
                Symbol(get(item, "kind", "dirichlet")),
                Float64(get(item, "value", 0.0)),
                Symbol(get(item, "method", "weak"))
            ))
    end
    return bcs
end

function _parse_initial_conditions(raw::AbstractDict)
    return InitialConditionConfig(
        Symbol(get(raw, "kind", "default")),
        _float_or_nothing(get(raw, "amplitude", nothing)),
        _float_or_nothing(get(raw, "value", nothing)),
        _float_or_nothing(get(raw, "mobile_amplitude", nothing)),
        _float_or_nothing(get(raw, "mobile_value", nothing)),
        _float_or_nothing(get(raw, "trap_occupancy", nothing))
    )
end

_float_or_nothing(value) = value === nothing ? nothing : Float64(value)

function _parse_parameters(cfg::AbstractDict)
    raw = get(cfg, "parameters", Dict{String, Any}())
    raw isa AbstractDict || throw(ArgumentError("parameters must be a TOML table"))

    params = Dict{Symbol, Float64}()
    for (key, value) in pairs(raw)
        if value isa Real
            params[Symbol(key)] = Float64(value)
        end
    end

    if !haskey(params, :t0)
        params[:t0] = 0.0
    end
    if !haskey(params, :tend)
        params[:tend] = 1.0
    end

    return (; pairs(params)...)
end

function _parse_saveat(value)
    value === nothing && return nothing
    value isa Real && return [Float64(value)]
    value isa AbstractVector ||
        throw(ArgumentError("solver.saveat must be a number, array, or omitted"))
    return Float64[Float64(v) for v in value]
end

function _parse_algorithm_spec(value)
    value isa AbstractString || throw(ArgumentError("solver.algorithm must be a string"))
    return Symbol(value)
end

function _parse_formulation_symbol(value)
    value isa AbstractString || throw(ArgumentError("solver.formulation must be a string"))
    return Symbol(lowercase(value))
end

function _normalize_problem_type(value)
    value == "toy_rd" && return "diffusion_1d"
    value == "toy_trapping" && return "trapping_1d"
    value == "fake_hotgates_trapping" && return "hotgates_trapping"
    return value
end

function _required_table(cfg::AbstractDict, key::AbstractString)
    haskey(cfg, key) || throw(ArgumentError("Missing required [$key] table"))
    value = cfg[key]
    value isa AbstractDict || throw(ArgumentError("[$key] must be a TOML table"))
    return value
end

function _required_string(cfg::AbstractDict, key::AbstractString)
    haskey(cfg, key) || throw(ArgumentError("Missing required key '$key'"))
    value = cfg[key]
    value isa AbstractString || throw(ArgumentError("'$key' must be a string"))
    return value
end
