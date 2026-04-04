function solution_dataframe(sol, model::SystemModel)
    layout = model.layout
    nx = model.context.nx
    x = model.context.mesh.x
    vars = variable_names(layout)

    rows = DataFrame(
        time = Float64[],
        x = Float64[],
        variable = String[],
        value = Float64[],
    )

    nvars = nvariables(layout)

    for (itime, t) in enumerate(sol.t)
        U = reshape(sol.u[itime], nvars, nx)
        for ivar in 1:nvars
            vname = String(vars[ivar])
            for ix in 1:nx
                push!(rows, (Float64(t), Float64(x[ix]), vname, Float64(U[ivar, ix])))
            end
        end
    end

    return rows
end

function compute_summary_timeseries(sol, model::SystemModel)
    layout = model.layout
    nx = model.context.nx
    dx = model.context.mesh.dx
    vars = variable_names(layout)
    nvars = nvariables(layout)

    df = DataFrame(time = Float64[])

    for v in vars
        df[!, Symbol("integral_" * String(v))] = Float64[]
        df[!, Symbol("max_" * String(v))] = Float64[]
        df[!, Symbol("min_" * String(v))] = Float64[]
    end

    for (itime, t) in enumerate(sol.t)
        U = reshape(sol.u[itime], nvars, nx)

        row = Dict{Symbol,Any}(:time => Float64(t))
        for ivar in 1:nvars
            vals = @view U[ivar, :]
            vname = String(vars[ivar])
            row[Symbol("integral_" * vname)] = sum(vals) * dx
            row[Symbol("max_" * vname)] = maximum(vals)
            row[Symbol("min_" * vname)] = minimum(vals)
        end
        push!(df, row)
    end

    return df
end

function solver_stats_dict(sol)
    stats = Dict{String,Any}()

    if hasproperty(sol, :retcode)
        stats["retcode"] = string(sol.retcode)
    end
    if hasproperty(sol, :stats)
        s = sol.stats
        for name in propertynames(s)
            stats[string(name)] = getproperty(s, name)
        end
    end

    return stats
end

function save_solution_hdf5(path::AbstractString, sol, model::SystemModel; run_info=nothing)
    layout = model.layout
    nx = model.context.nx
    x = model.context.mesh.x
    vars = string.(variable_names(layout))
    nvars = nvariables(layout)
    nt = length(sol.t)

    data = Array{Float64}(undef, nvars, nx, nt)
    for it in 1:nt
        data[:, :, it] .= reshape(sol.u[it], nvars, nx)
    end

    h5open(path, "w") do h5
        h5["/time"] = collect(Float64.(sol.t))
        h5["/mesh/x"] = collect(Float64.(x))
        h5["/fields/data"] = data
        h5["/fields/variable_names"] = vars

        attrs(h5)["nx"] = nx
        attrs(h5)["nvars"] = nvars

        if run_info !== nothing
            attrs(h5)["run_id"] = run_info.run_id
            attrs(h5)["start_time"] = string(run_info.start_time)
            if run_info.config_path !== nothing
                attrs(h5)["config_path"] = run_info.config_path
            end
        end

        stats = solver_stats_dict(sol)
        for (k, v) in stats
            attrs(h5)["solver_" * k] = string(v)
        end
    end

    return path
end

function save_summary_csv(path::AbstractString, sol, model::SystemModel)
    df = compute_summary_timeseries(sol, model)
    CSV.write(path, df)
    return path
end

function save_diagnostics_csv(path::AbstractString, sol)
    df = DataFrame(
        time = Float64.(sol.t),
    )
    CSV.write(path, df)
    return path
end
