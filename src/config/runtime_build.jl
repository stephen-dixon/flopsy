"""
    build_simulation(ctx, problem_name = _default_problem_name(ctx))

Assemble a named problem from a validated `BuildContext`.
"""
function build_simulation(ctx::BuildContext, problem_name::Symbol = _default_problem_name(ctx))
    haskey(ctx.problems, problem_name) ||
        throw(ConfigValidationError("Unknown problem block '$problem_name'"))

    problem_def = ctx.problems[problem_name]
    mesh = _lookup_named(ctx.meshes, problem_def.mesh, :mesh)
    backend = _lookup_named(ctx.backends, problem_def.backend, :backend)
    ics = ICDefinition[_lookup_named(ctx.ics, name, :ic) for name in problem_def.ics]
    bcs = BCDefinition[_lookup_named(ctx.bcs, name, :bc) for name in problem_def.bcs]
    outputs = OutputDefinition[_lookup_named(ctx.outputs, name, :output)
                               for name in problem_def.outputs]

    temperature = if problem_def.temperature !== nothing
        _lookup_named(ctx.temperatures, problem_def.temperature, :temperature)
    else
        nothing
    end

    _validate_bc_targets!(backend, bcs)
    _validate_ic_targets!(backend, ics)

    model = backend.build_model(mesh, bcs, temperature)
    u0 = zeros(Float64, nvariables(model.layout) * model.context.nx)
    for ic in ics
        ic.apply!(u0, backend, model)
    end

    solver_config = _build_solver_config(problem_def.solver)
    sim_problem = SimulationProblem(
        model, u0, problem_def.tspan, solver_config, problem_def)
    return ConfiguredSimulation(sim_problem, outputs, ctx)
end

function solve(plan::ConfiguredSimulation; output_dir = nothing)
    return _solve_configured_simulation(plan; output_dir = output_dir)
end

function _solve_configured_simulation(plan::ConfiguredSimulation; output_dir = nothing)
    result = solve(plan.problem)
    for output in plan.outputs
        _materialize_output!(result, output; output_dir = output_dir)
    end
    return result
end

"""
    validate_input_deck(path; registry = build_registry(), problem_name = nothing)

Parse, validate, and assemble an input deck without running the solver.
"""
function validate_input_deck(
        path::AbstractString;
        registry::SyntaxRegistry = build_registry(),
        problem_name = nothing
)
    deck = parse_input_deck(path)
    ctx = build_context(deck; registry = registry)
    pname = problem_name === nothing ? _default_problem_name(ctx) : Symbol(problem_name)
    build_simulation(ctx, pname)
    return ctx
end

"""
    run_input_deck(path; registry = build_registry(), problem_name = nothing, output_dir = nothing)

Parse, validate, assemble, and run a registry-driven input deck.

- `problem_name` — select a specific `[problem.<name>]` block (default: alphabetically first)
- `output_dir`   — override output paths to write all outputs under this directory
"""
function run_input_deck(
        path::AbstractString;
        registry::SyntaxRegistry = build_registry(),
        problem_name = nothing,
        output_dir = nothing
)
    deck = parse_input_deck(path)
    ctx = build_context(deck; registry = registry)
    pname = problem_name === nothing ? _default_problem_name(ctx) : Symbol(problem_name)
    return solve(build_simulation(ctx, pname); output_dir = output_dir)
end

function _default_problem_name(ctx::BuildContext)
    isempty(ctx.problems) &&
        throw(ConfigValidationError("No [problem.<name>] blocks were defined"))
    return only(sort!(collect(keys(ctx.problems)); by = string))
end

function _lookup_named(dict::Dict{Symbol, Any}, name::Symbol, domain::Symbol)
    haskey(dict, name) ||
        throw(ConfigValidationError("Unknown $(domain) reference '$name'"))
    return dict[name]
end

function _validate_bc_targets!(backend, bcs::AbstractVector)
    for bc in bcs
        info = _species_info(backend, bc.species)
        info.boundary_target ||
            throw(ConfigValidationError(
                "Boundary condition `$(bc.name)` targets species `$(bc.species)`, which does not support boundary treatment",
            ))
        info.transport == :diffusive ||
            throw(ConfigValidationError(
                "Boundary condition `$(bc.name)` targets species `$(bc.species)`, which is not diffusive",
            ))
    end
end

function _validate_ic_targets!(backend, ics)
    assigned = Set{Symbol}()
    for ic in ics
        for name in ic.affects
            _species_info(backend, name)
            name in assigned &&
                throw(ConfigValidationError("Overlapping initial conditions detected for species `$(name)`"))
            push!(assigned, name)
        end
    end
end

function _species_info(backend, name::Symbol)
    for info in backend.species
        info.name == name && return info
    end
    available = join(string.(getfield.(backend.species, :name)), ", ")
    throw(ConfigValidationError("Unknown species `$(name)`. Available species: $(available)"))
end

_species_index(backend, name::Symbol) = _species_info(backend, name).index

function _materialize_output!(result::SimulationResult, output::OutputDefinition; output_dir = nothing)
    file = _resolve_output_path(output.file, output_dir)

    if output.type_name == :hdf5
        write_field_output_hdf5(result, file;
            export_equilibrium_trapped = output.export_equilibrium_trapped,
            export_retention_total = output.export_retention_total,
            export_retention_by_occupation = output.export_retention_by_occupation,
            export_Deff = output.export_Deff)
        if output.xdmf_path !== nothing
            xdmf = _resolve_output_path(output.xdmf_path, output_dir)
            write_xdmf_for_flopsy_h5(file, xdmf)
        end
        return file
    elseif output.type_name == :summary_csv
        write_summary_csv(result, file)
        return file
    end
    throw(ConfigValidationError("Unsupported output type $(output.type_name)"))
end

function _resolve_output_path(path::AbstractString, output_dir)
    output_dir === nothing && return String(path)
    return joinpath(String(output_dir), basename(String(path)))
end
