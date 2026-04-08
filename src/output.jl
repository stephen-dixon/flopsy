using CSV
using DataFrames
using HDF5

function wrap_result(model::SystemModel, sol, config;
    summaries::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    metadata::Dict{String,Any}=Dict{String,Any}(),
)
    return SimulationResult(model, sol, config, summaries, metadata)
end


"""
    variable_snapshot(result, var; time_index=length(result.solution.u))

Return the spatial profile of variable `var` at a given saved time index.

`var` may be:
- an integer variable index
- a Symbol variable name
"""
function variable_snapshot(result::SimulationResult, var; time_index=length(result.solution.u))
    model = result.model
    layout = model.layout
    nx = model.context.nx

    ivar = _resolve_variable_index(layout, var)
    U = state_view(result.solution.u[time_index], layout, nx)
    return copy(@view U[ivar, :])
end


"""
    variable_timeseries(result, var)

Return a matrix of size (nt, nx) for one variable over all saved times.
"""
function variable_timeseries(result::SimulationResult, var)
    model = result.model
    layout = model.layout
    nx = model.context.nx
    ivar = _resolve_variable_index(layout, var)

    nt = length(result.solution.u)
    out = zeros(eltype(result.solution.u[1]), nt, nx)

    for it in 1:nt
        U = state_view(result.solution.u[it], layout, nx)
        out[it, :] .= @view U[ivar, :]
    end

    return out
end


"""
    integrated_variable(result, var)

Compute a simple spatial integral over time for variable `var`.
Returns a vector with one value per saved time.
"""
function integrated_variable(result::SimulationResult, var)
    model = result.model
    mesh = model.context.mesh
    vals = variable_timeseries(result, var)
    nt = size(vals, 1)

    out = zeros(eltype(vals), nt)
    for it in 1:nt
        out[it] = sum(@view vals[it, :]) * mesh.dx
    end
    return out
end


function build_summary_dataframe(result::SimulationResult)
    layout = result.model.layout
    names = variable_names(layout)
    t = result.solution.t

    df = DataFrame(time = t)

    for name in names
        df[!, Symbol("integral_", name)] = integrated_variable(result, name)
    end

    if haskey(result.summaries, :extra_timeseries)
        extra = result.summaries[:extra_timeseries]
        for (k, v) in pairs(extra)
            df[!, Symbol(k)] = v
        end
    end

    return df
end


"""
    write_summary_csv(result, path)

Write summary/integrated quantities over time to CSV.
"""
function write_summary_csv(result::SimulationResult, path::AbstractString)
    df = build_summary_dataframe(result)
    CSV.write(path, df)
    return path
end


"""
    write_field_output_hdf5(result, path)

Write pointwise variable outputs over time to HDF5.

Layout:
- /time                  : saved times
- /mesh/x                : spatial coordinates
- /fields/<varname>      : matrix (nt, nx)
- /metadata/...          : selected metadata as string attributes/datasets
"""
function write_field_output_hdf5(result::SimulationResult, path::AbstractString)
    model = result.model
    layout = model.layout
    mesh = model.context.mesh
    t = result.solution.t

    h5open(path, "w") do h5
        h5["time"] = collect(t)
        h5["mesh/x"] = mesh.x
        h5["mesh/dx"] = mesh.dx

        g_fields = create_group(h5, "fields")
        for name in variable_names(layout)
            g_fields[string(name)] = variable_timeseries(result, name)
        end

        g_meta = create_group(h5, "metadata")
        for (k, v) in pairs(result.metadata)
            try
                g_meta[string(k)] = string(v)
            catch
                # fallback silently for now
            end
        end
    end

    return path
end


function print_run_banner(config, solver_config::SolverConfig, model::SystemModel)
    println("============================================================")
    println("Flopsy simulation start")
    println("------------------------------------------------------------")
    versions = library_versions()
    println("Julia version       : ", versions["julia_version"])
    println("Flopsy version      : ", versions["flopsy_version"])
    println("------------------------------------------------------------")
    println("Variables           : ", join(string.(variable_names(model.layout)), ", "))
    println("Number of variables : ", nvariables(model.layout))
    println("Number of nodes     : ", model.context.nx)
    println("Domain length       : ", model.context.mesh.x[end] - model.context.mesh.x[1])
    println("dx                  : ", model.context.mesh.dx)
    println("Formulation         : ", typeof(solver_config.formulation))
    println("Algorithm           : ", solver_config.algorithm)
    println("AbsTol              : ", solver_config.abstol)
    println("RelTol              : ", solver_config.reltol)
    println("============================================================")
end


function _resolve_variable_index(layout::VariableLayout, var::Integer)
    1 <= var <= nvariables(layout) || throw(ArgumentError("Variable index $var out of range"))
    return var
end

function _resolve_variable_index(layout::VariableLayout, var::Symbol)
    names = variable_names(layout)
    idx = findfirst(==(var), names)
    idx === nothing && throw(ArgumentError("Unknown variable name $var"))
    return idx
end

# optional plotting helpers if Plots.jl is added later
# using Plots

# function plot_snapshot(result::SimulationResult, var; time_index=length(result.solution.u))
#     x = result.model.context.mesh.x
#     y = variable_snapshot(result, var; time_index=time_index)
#     return plot(x, y, xlabel="x", ylabel=string(var), label=string(var))
# end

function library_versions()
    return Dict(
        "julia_version" => string(VERSION),
        "flopsy_version" => "0.1.0-dev",
        "scimlbase_loaded" => isdefined(Main, :SciMLBase) || isdefined(Flopsy, :SciMLBase),
        "ordinarydiffeq_loaded" => isdefined(Main, :OrdinaryDiffEq) || isdefined(Flopsy, :OrdinaryDiffEq),
        "sundials_loaded" => isdefined(Main, :Sundials) || isdefined(Flopsy, :Sundials),
    )
end

function solver_stats_dict(result::SimulationResult)
    sol = result.solution
    stats = Dict{String,Any}()

    for name in propertynames(sol.stats)
        try
            stats[string(name)] = getproperty(sol.stats, name)
        catch
        end
    end

    return stats
end
