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
        summaries::Dict{Symbol, Any} = Dict{Symbol, Any}(),
        metadata::Dict{String, Any} = Dict{String, Any}()
)
    return SimulationResult(model, sol, config, summaries, metadata)
end

"""
    variable_snapshot(result, var; time_index=length(result.solution.u))

Return the spatial profile of variable `var` at a given saved time index.

`var` may be an integer variable index or a Symbol variable name.
"""
function variable_snapshot(result::SimulationResult, var; time_index = length(result.solution.u))
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
        df[!, Symbol("left_flux_", varname)] = fl.left
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
    write_xdmf_for_flopsy_h5(h5_path, xdmf_path=nothing) -> String

Generate an XDMF companion file for a Flopsy HDF5 field output file so the
result can be opened directly in ParaView.
"""
function write_xdmf_for_flopsy_h5(h5_path::AbstractString, xdmf_path = nothing)
    out = xdmf_path === nothing ? replace(h5_path, r"\.h5$" => ".xdmf") : String(xdmf_path)
    h5_ref = basename(h5_path)

    times = Float64[]
    x = Float64[]
    field_names = String[]
    h5open(h5_path, "r") do h5
        times = read(h5["time"])
        x = read(h5["mesh/x"])
        field_names = sort!(collect(String.(keys(h5["fields"]))))
    end

    nt = length(times)
    nx = length(x)

    open(out, "w") do io
        println(io, "<?xml version=\"1.0\" ?>")
        println(io, "<Xdmf Version=\"3.0\">")
        println(io, "  <Domain>")
        println(io, "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">")
        for (it, tval) in enumerate(times)
            println(io, "      <Grid Name=\"step_$(it - 1)\" GridType=\"Uniform\">")
            println(io, "        <Time Value=\"$(Float64(tval))\"/>")
            println(io, "        <Topology TopologyType=\"1DSMesh\" NumberOfElements=\"$(nx)\"/>")
            println(io, "        <Geometry GeometryType=\"X\">")
            println(io,
                "          <DataItem Format=\"HDF\" Dimensions=\"$(nx)\" NumberType=\"Float\" Precision=\"8\">$(h5_ref):/mesh/x</DataItem>")
            println(io, "        </Geometry>")
            for name in field_names
                println(io, "        <Attribute Name=\"$(name)\" AttributeType=\"Scalar\" Center=\"Node\">")
                println(io,
                    "          <DataItem ItemType=\"HyperSlab\" Dimensions=\"$(nx)\" Type=\"HyperSlab\">")
                println(io, "            <DataItem Dimensions=\"3 2\" Format=\"XML\">")
                println(io, "              $(it - 1) 0")
                println(io, "              1 1")
                println(io, "              1 $(nx)")
                println(io, "            </DataItem>")
                println(io,
                    "            <DataItem Format=\"HDF\" Dimensions=\"$(nt) $(nx)\" NumberType=\"Float\" Precision=\"8\">$(h5_ref):/fields/$(name)</DataItem>")
                println(io, "          </DataItem>")
                println(io, "        </Attribute>")
            end
            println(io, "      </Grid>")
        end
        println(io, "    </Grid>")
        println(io, "  </Domain>")
        println(io, "</Xdmf>")
    end

    return out
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
        time_index::Union{Int, Symbol} = :last
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
    check_mass_conservation(result; rtol=1e-3) -> NamedTuple

Verify that total hydrogen is conserved up to the integrated surface flux.

Computes the total hydrogen inventory at each saved time:
    H(t) = ∑ᵢ ∫₀ᴸ uᵢ(x,t) dx    (sum over all variables, spatial integral)

and the cumulative surface loss via the trapezoidal rule:
    Φ(t) = ∫₀ᵗ [left_flux(s) + right_flux(s)] ds

Mass balance: H(t) + Φ(t) ≈ H(0) for all saved times.

The surface flux is computed using the mass-conserving formula:
    left_flux  = D * (U[1]  - g_left(t))  / dx
    right_flux = D * (U[nx] - g_right(t)) / dx
where g(t) is the Dirichlet BC value.  For closed systems (Neumann / no boundary
operator), Φ(t) = 0 and conservation reduces to checking H is constant.

# Returns
A `NamedTuple` with fields:
- `conserved`            — `true` if the maximum relative balance error ≤ `rtol`
- `max_relative_error`   — maximum relative error over all saved times
- `total_hydrogen`       — `Vector{Float64}` of H(t)
- `cumulative_flux`      — `Vector{Float64}` of Φ(t)
- `balance`              — `H(t) + Φ(t)` (should be constant ≈ H(0))

# Notes
- Trapping models satisfy global mass conservation exactly in continuous time;
  the numerical residual reflects time-integration and spatial-quadrature errors.
- Only `WeakDirichletBoundaryOperator` and `OperatorSum` thereof contribute surface
  flux.  Penalty / mass-matrix / callback / eliminated methods are not yet accounted
  for (these are rarely combined with mass-conservation diagnostics).
"""
function check_mass_conservation(result::SimulationResult; rtol::Real = 1e-3)
    layout = result.model.layout
    t = result.solution.t
    nt = length(t)

    # Total hydrogen at each saved time
    H = zeros(Float64, nt)
    for name in variable_names(layout)
        H .+= integrated_variable(result, name)
    end

    # Compute mass-conserving surface flux (correct Dirichlet formula D*(U[bdy]-g)/dx)
    fluxes = _mass_conserving_boundary_fluxes(result)
    total_flux = zeros(Float64, nt)
    for (_, fl) in pairs(fluxes)
        total_flux .+= fl.left .+ fl.right
    end

    # Cumulative integral via trapezoidal rule
    cum_flux = zeros(Float64, nt)
    for k in 2:nt
        dt = t[k] - t[k - 1]
        cum_flux[k] = cum_flux[k - 1] + 0.5 * dt * (total_flux[k] + total_flux[k - 1])
    end

    H0 = H[1]
    balance = H .+ cum_flux          # should equal H0 at all times
    denom = max(abs(H0), 1e-30)
    balance_errors = abs.(balance .- H0) ./ denom
    max_rel_error = maximum(balance_errors)

    return (
        conserved = max_rel_error ≤ rtol,
        max_relative_error = max_rel_error,
        total_hydrogen = H,
        cumulative_flux = cum_flux,
        balance = balance
    )
end

"""
    _mass_conserving_boundary_fluxes(result) -> Dict{Symbol, NamedTuple}

Compute boundary fluxes using the mass-conserving formula `D*(U[boundary]-g)/dx`.
Only `WeakDirichletBoundaryOperator` contributes; Neumann (closed) systems return
an empty dict (zero flux, H conserved trivially).

This differs from `surface_diffusive_fluxes`, which uses `D*(U[2]-U[1])/dx`
(interior gradient approximation, suitable for TDS flux plots but not for mass
balance checks).
"""
function _mass_conserving_boundary_fluxes(result::SimulationResult)
    model = result.model
    layout = model.layout
    nx = model.context.nx
    dx = model.context.mesh.dx
    ctx = model.context
    t_arr = result.solution.t
    nt = length(t_arr)
    names = variable_names(layout)

    fluxes = Dict{Symbol, NamedTuple}()

    bcop = model.operators.boundary
    bcop === nothing && return fluxes

    _accumulate_weak_dirichlet_fluxes!(
        fluxes, bcop, result, layout, nx, dx, ctx, t_arr, nt, names)
    return fluxes
end

function _accumulate_weak_dirichlet_fluxes!(fluxes, op::WeakDirichletBoundaryOperator,
        result, layout, nx, dx, ctx, t_arr, nt, names)
    # Neither side has a BC → nothing to add
    op.left === nothing && op.right === nothing && return

    vars = op.selector(layout)
    for ivar in vars
        fl = get!(fluxes, names[ivar]) do
            (left = zeros(Float64, nt), right = zeros(Float64, nt))
        end

        for it in 1:nt
            t = t_arr[it]
            T_val = op.temperature !== nothing ?
                    Float64(temperature_at(op.temperature, ctx, t, 1)) : NaN
            D = Flopsy._eval_D(op.coefficients, ivar, 1, T_val)
            U = state_view(result.solution.u[it], layout, nx)

            if op.left !== nothing
                fl.left[it] += D * (U[ivar, 1] - op.left(t)) / dx
            end
            if op.right !== nothing
                fl.right[it] += D * (U[ivar, nx] - op.right(t)) / dx
            end
        end
    end
end

# OperatorSum as boundary: recurse into sub-operators
function _accumulate_weak_dirichlet_fluxes!(fluxes, op::OperatorSum,
        result, layout, nx, dx, ctx, t_arr, nt, names)
    for sub in op.ops
        _accumulate_weak_dirichlet_fluxes!(
            fluxes, sub, result, layout, nx, dx, ctx, t_arr, nt, names)
    end
end

# Default: ignore operators that don't contribute a mass-conserving boundary flux
function _accumulate_weak_dirichlet_fluxes!(fluxes, op::AbstractOperator,
        result, layout, nx, dx, ctx, t_arr, nt, names)
    return
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
    1 <= var <= nvariables(layout) ||
        throw(ArgumentError("Variable index $var out of range"))
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
        "julia_version" => string(VERSION),
        "flopsy_version" => string(pkgversion(@__MODULE__)),
        "scimlbase_loaded" => isdefined(Main, :SciMLBase) || isdefined(Flopsy, :SciMLBase),
        "ordinarydiffeq_loaded" => isdefined(Main, :OrdinaryDiffEq) ||
                                   isdefined(Flopsy, :OrdinaryDiffEq),
        "sundials_loaded" => isdefined(Main, :Sundials) || isdefined(Flopsy, :Sundials)
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
    hasproperty(sol, :stats) && sol.stats !== nothing || return Dict{String, Any}()

    stats = Dict{String, Any}()
    for name in propertynames(sol.stats)
        try
            stats[string(name)] = getproperty(sol.stats, name)
        catch
        end
    end

    return stats
end
