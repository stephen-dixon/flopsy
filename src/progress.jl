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
