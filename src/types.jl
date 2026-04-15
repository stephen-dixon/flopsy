"Abstract base type for all system models."
abstract type AbstractSystemModel end

"Abstract base type for all operators (reaction, diffusion, constraint)."
abstract type AbstractOperator end

"Abstract type for reaction operators contributing to `du/dt`."
abstract type AbstractReactionOperator <: AbstractOperator end

"Abstract type for diffusion operators contributing spatial transport terms."
abstract type AbstractDiffusionOperator <: AbstractOperator end

"Abstract type for algebraic constraint operators (DAE path)."
abstract type AbstractConstraintOperator <: AbstractOperator end

"Abstract type for time-integration formulations."
abstract type AbstractFormulation end

"""
    Mesh1D(length_domain, nx)

Uniform 1D mesh with `nx` nodes spanning `[0, length_domain]`.

Fields:
- `x`  — node positions (`Vector{Float64}`, length `nx`)
- `dx` — uniform node spacing
"""
struct Mesh1D{T}
    x::Vector{T}
    dx::T
end

struct RunInfo
    run_id::String
    start_time::DateTime
    config_path::Union{Nothing, String}
    package_versions::Dict{String, String}
end

"""
    SimulationProblem

Typed container representing a fully assembled simulation ready to execute.
This is the hand-off boundary between the config/problem-template layer and the
core solver framework.
"""
struct SimulationProblem{M, U, T, S, C}
    model::M
    u0::U
    tspan::T
    solver_config::S
    config::C
end

"""
    SystemContext

Runtime context passed to every operator call.  Holds the mesh, node count,
variable layout, auxiliary data, and scratch buffers.
"""
struct SystemContext{L, M, A, S}
    layout::L
    nx::Int
    mesh::M
    aux::A
    scratch::S
end

"""
    SystemModel

Assembled model combining a `VariableLayout`, operator NamedTuple, and `SystemContext`.
"""
struct SystemModel{L, O, C} <: AbstractSystemModel
    layout::L
    operators::O
    context::C
end

"""
    SolverConfig(; formulation, algorithm, abstol, reltol, saveat, ...)

Configuration for the time integrator.

- `formulation` — e.g. `UnsplitFormulation()`
- `algorithm`   — SciML algorithm, e.g. `Rodas5(autodiff=AutoFiniteDiff())`
- `abstol`, `reltol` — solver tolerances
- `saveat`      — times to save solution (nothing = solver default)
- `dt`          — initial step size hint (nothing = auto)
- `kwargs`      — extra keyword arguments forwarded to `SciMLBase.solve`
"""
Base.@kwdef struct SolverConfig{F, A, AT, RT, S, D, K}
    formulation::F
    algorithm::A
    abstol::AT = 1e-8
    reltol::RT = 1e-6
    saveat::S = nothing
    dt::D = nothing
    show_progress::Bool = true
    show_solver_stats::Bool = true
    write_convergence_trace::Bool = false
    kwargs::K = NamedTuple()
end

"""
    SimulationResult

Bundles the model, SciML solution, configuration, and optional summary/metadata
from a completed simulation.  Access fields with the output helper functions:
`variable_timeseries`, `variable_snapshot`, `integrated_variable`,
`surface_diffusive_fluxes`, `build_summary_dataframe`, etc.
"""
struct SimulationResult{M, S, C}
    model::M
    solution::S
    config::C
    summaries::Dict{Symbol, Any}
    metadata::Dict{String, Any}
end
