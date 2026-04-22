function package_versions()
    return Dict(
        "Flopsy" => "0.1.0-dev",
        "OrdinaryDiffEq" => string(Base.PkgId(OrdinaryDiffEq).version),
        "SciMLBase" => string(Base.PkgId(SciMLBase).version),
        "HDF5" => string(Base.PkgId(HDF5).version),
        "CSV" => string(Base.PkgId(CSV).version),
        "DataFrames" => string(Base.PkgId(DataFrames).version)
    )
end

function make_run_info(; config_path::Union{Nothing, String} = nothing)
    return RunInfo(
        string(uuid4()),
        now(),
        config_path,
        package_versions()
    )
end

function print_run_banner(
        model::SystemModel, u0, tspan, solver_config::SolverConfig, run_info::RunInfo)
    println("==================================================")
    println("Flopsy run starting")
    println("Run ID:         ", run_info.run_id)
    println("Start time:     ", run_info.start_time)
    println("Config path:    ", something(run_info.config_path, "<none>"))
    println("Formulation:    ", typeof(solver_config.formulation))
    println("Algorithm:      ", typeof(solver_config.algorithm))
    println("tspan:          ", tspan)
    println("nx:             ", model.context.nx)
    println("nvars:          ", nvariables(model.layout))
    println("ndofs:          ", length(u0))
    println("Variable names: ", join(string.(variable_names(model.layout)), ", "))
    println("Package versions:")
    for (k, v) in sort(collect(run_info.package_versions); by = first)
        println("  ", rpad(k, 18), v)
    end
    println("==================================================")
end

mutable struct ProgressReporter
    enabled::Bool
    io::IO
    label::String
    t0::Float64
    t_end::Float64
    width::Int
    last_percent::Int
    last_print_ns::UInt64
    started_ns::UInt64
    finished::Bool
end

function _progress_enabled(show_progress::Bool)
    show_progress || return false
    mode = lowercase(get(ENV, "FLOPSY_PROGRESS", "auto"))
    mode in ("0", "false", "off", "no") && return false
    mode in ("1", "true", "on", "yes", "force") && return true
    return stderr isa Base.TTY
end

function _start_progress(tspan; show_progress::Bool = true, label::AbstractString = "Solving")
    t0 = Float64(tspan[1])
    t_end = Float64(tspan[2])
    now = time_ns()
    reporter = ProgressReporter(
        _progress_enabled(show_progress),
        stderr,
        String(label),
        t0,
        t_end,
        32,
        -1,
        now,
        now,
        false
    )
    _update_progress!(reporter, t0; force = true)
    return reporter
end

function _update_progress!(p::ProgressReporter, t; force::Bool = false)
    p.enabled || return nothing
    p.finished && return nothing

    denom = p.t_end - p.t0
    frac = denom == 0 ? 1.0 : clamp((Float64(t) - p.t0) / denom, 0.0, 1.0)
    percent = Int(floor(100 * frac))
    now = time_ns()
    if !force && percent == p.last_percent && now - p.last_print_ns < UInt64(200_000_000)
        return nothing
    end

    filled = Int(floor(frac * p.width))
    empty = p.width - filled
    elapsed = (now - p.started_ns) / 1.0e9
    print(p.io,
        "\r",
        p.label,
        " [",
        repeat("=", filled),
        repeat(" ", empty),
        "] ",
        lpad(string(percent), 3),
        "%  t=",
        _format_progress_time(Float64(t)),
        "/",
        _format_progress_time(p.t_end),
        "  elapsed=",
        _format_elapsed(elapsed)
    )
    flush(p.io)
    p.last_percent = percent
    p.last_print_ns = now
    return nothing
end

function _finish_progress!(p::ProgressReporter)
    p.enabled || return nothing
    p.finished && return nothing
    p.last_percent < 100 && _update_progress!(p, p.t_end; force = true)
    println(p.io)
    p.finished = true
    return nothing
end

function _progress_callback(p::ProgressReporter)
    p.enabled || return nothing
    condition(u, t, integrator) = true
    affect!(integrator) = _update_progress!(p, integrator.t)
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

function _merge_callback(kwargs::NamedTuple, cb)
    cb === nothing && return kwargs
    existing = haskey(kwargs, :callback) ? kwargs.callback : nothing
    callback = existing === nothing ? cb : CallbackSet(existing, cb)
    return merge(kwargs, (; callback = callback))
end

function _format_progress_time(t::Float64)
    if isinteger(t)
        return string(Int(t))
    elseif abs(t) >= 1e4 || (abs(t) < 1e-3 && t != 0)
        return string(round(t; sigdigits = 4))
    else
        return string(round(t; digits = 3))
    end
end

function _format_elapsed(seconds::Real)
    seconds < 60 && return string(round(seconds; digits = 1), "s")
    minutes = floor(Int, seconds / 60)
    rem = floor(Int, seconds - 60 * minutes)
    return string(minutes, "m", lpad(string(rem), 2, "0"), "s")
end
