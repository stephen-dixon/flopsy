using CSV
using DataFrames
using HDF5

"""
    wrap_result(model, sol, config; summaries=Dict(), metadata=Dict()) -> SimulationResult

Wrap a SciML solution object together with the model, configuration, and optional
summary/metadata dicts into a `SimulationResult`.

`summaries` may contain a `:extra_timeseries` key mapping string names to
`Vector{Float64}` — these appear as additional columns in the summary CSV.
"""
function wrap_result(model::SystemModel, sol, config;
    summaries::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    metadata::Dict{String,Any}=Dict{String,Any}(),
)
    return SimulationResult(model, sol, config, summaries, metadata)
end


"""
    variable_snapshot(result, var; time_index=length(result.solution.u))

Return the spatial profile of variable `var` at a given saved time index.

`var` may be an integer variable index or a Symbol variable name.
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

Return a matrix of size `(nt, nx)` for one variable over all saved times.
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

Spatial integral of variable `var` over time (simple trapezoidal rule with uniform dx).
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


"""
    surface_diffusive_fluxes(result) -> Dict{Symbol, NamedTuple}

Compute the outward diffusive flux at both domain boundaries for each diffusing
variable over all saved times.

Returns a `Dict` keyed by variable name.  Each entry is a `NamedTuple`:
    `(left = Vector{Float64}, right = Vector{Float64})`

Sign convention: **positive in the direction of the outward-facing normal**.

- **left**:  `D * (C[2] - C[1]) / dx`     — outward normal is −x; equals `D * ∂C/∂x|_{x=0}`
- **right**: `D * (C[nx-1] - C[nx]) / dx` — outward normal is +x; equals `−D * ∂C/∂x|_{x=L}`

Both are positive during TDS desorption where the surface concentration is lower than
the bulk (vacuum boundary conditions), giving a positive outward flux at both faces.

Returns an empty `Dict` if the model has no diffusion operator, or if the operator
type does not implement `surface_fluxes`.
"""
function surface_diffusive_fluxes(result::SimulationResult)
    diffop = result.model.operators.diffusion
    diffop === nothing && return Dict{Symbol, NamedTuple}()
    return surface_fluxes(diffop, result)
end


"""
    build_summary_dataframe(result)

Build a `DataFrame` of scalar time-series summaries:
- `time`                      — saved times
- `integral_<var>`            — spatial integral of each variable
- `left_flux_<var>`           — left-surface desorption flux for each diffusing variable
- `right_flux_<var>`          — right-surface desorption flux for each diffusing variable
- any extra entries from `result.summaries[:extra_timeseries]`
"""
function build_summary_dataframe(result::SimulationResult)
    layout = result.model.layout
    names = variable_names(layout)
    t = result.solution.t

    df = DataFrame(time = t)

    for name in names
        df[!, Symbol("integral_", name)] = integrated_variable(result, name)
    end

    fluxes = surface_diffusive_fluxes(result)
    for (varname, fl) in pairs(fluxes)
        df[!, Symbol("left_flux_",  varname)] = fl.left
        df[!, Symbol("right_flux_", varname)] = fl.right
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
- `/time`              — saved times (Float64 vector)
- `/mesh/x`            — spatial node coordinates
- `/mesh/dx`           — node spacing
- `/fields/<varname>`  — matrix `(nt, nx)` for each variable
- `/metadata/...`      — selected metadata as string datasets
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
            end
        end
    end

    return path
end


"""
    load_ic_from_hdf5(path, model; time_index=:last) -> Vector{Float64}

Load an initial condition from a previously written HDF5 field output file.

Matches variables by name between the file and `model.layout`.  Variables present
in the model but absent from the file are initialised to zero with a warning.
This allows chaining simulations where the variable set differs (e.g. adding
a new trap type) or the number of nodes changes (mesh must match).

# Arguments
- `path`        — path to the HDF5 file written by `write_field_output_hdf5`
- `model`       — the target `SystemModel`
- `time_index`  — integer index into the saved time axis, or `:last` (default)
                  to use the final saved state

# Example — chaining implantation and desorption runs
```julia
# Implantation run
res_implant = run_simulation("implant.toml")
write_field_output_hdf5(res_implant, "implant_out.h5")

# Desorption run starting from the end of implantation
u0 = load_ic_from_hdf5("implant_out.h5", model_desorption)
sol = solve_problem(model_desorption, u0, (0.0, t_tds), solver_config)
```
"""
function load_ic_from_hdf5(
    path::AbstractString,
    model::SystemModel;
    time_index::Union{Int,Symbol} = :last,
)
    layout = model.layout
    nx = model.context.nx
    u0 = zeros(Float64, nvariables(layout) * nx)
    U0 = state_view(u0, layout, nx)

    h5open(path, "r") do h5
        available = Set(keys(h5["fields"]))

        for (ivar, name) in enumerate(variable_names(layout))
            sname = string(name)

            if sname ∈ available
                data = read(h5["fields"][sname])   # (nt, nx)
                size(data, 2) == nx || throw(DimensionMismatch(
                    "Variable $name in $path has $(size(data,2)) nodes; model expects $nx"
                ))
                idx = time_index === :last ? size(data, 1) : Int(time_index)
                U0[ivar, :] .= data[idx, :]
            else
                @warn "Variable $name not found in $path — initialising to zero"
            end
        end
    end

    return u0
end


"""
    print_run_banner(config, solver_config, model)

Print a human-readable summary of the simulation setup to stdout.
"""
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


function library_versions()
    return Dict(
        "julia_version"         => string(VERSION),
        "flopsy_version"        => string(pkgversion(@__MODULE__)),
        "scimlbase_loaded"      => isdefined(Main, :SciMLBase) || isdefined(Flopsy, :SciMLBase),
        "ordinarydiffeq_loaded" => isdefined(Main, :OrdinaryDiffEq) || isdefined(Flopsy, :OrdinaryDiffEq),
        "sundials_loaded"       => isdefined(Main, :Sundials) || isdefined(Flopsy, :Sundials),
    )
end

"""
    solver_stats_dict(result) -> Dict{String,Any}

Extract solver statistics from the SciML solution into a plain dictionary
(number of steps, function evaluations, etc.).
"""
function solver_stats_dict(result::SimulationResult)
    sol = result.solution

    # SplitSolution and other custom solution types may not have a .stats field.
    hasproperty(sol, :stats) && sol.stats !== nothing || return Dict{String,Any}()

    stats = Dict{String,Any}()
    for name in propertynames(sol.stats)
        try
            stats[string(name)] = getproperty(sol.stats, name)
        catch
        end
    end

    return stats
end
