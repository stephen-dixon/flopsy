abstract type AbstractSystemModel end

abstract type AbstractOperator end
abstract type AbstractReactionOperator <: AbstractOperator end
abstract type AbstractDiffusionOperator <: AbstractOperator end
abstract type AbstractConstraintOperator <: AbstractOperator end

abstract type AbstractFormulation end

struct Mesh1D{T}
    x::Vector{T}
    dx::T
end

struct RunInfo
    run_id::String
    start_time::DateTime
    config_path::Union{Nothing,String}
    package_versions::Dict{String,String}
end

struct SystemContext{L,M,A,S}
    layout::L
    nx::Int
    mesh::M
    aux::A
    scratch::S
end

struct SystemModel{L,O,C} <: AbstractSystemModel
    layout::L
    operators::O
    context::C
end

Base.@kwdef struct SolverConfig
    formulation::AbstractFormulation
    algorithm
    abstol::Float64 = 1e-8
    reltol::Float64 = 1e-6
    saveat = nothing
    dt = nothing
    show_progress::Bool = true
    show_solver_stats::Bool = true
    write_convergence_trace::Bool = false
    kwargs::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

struct SimulationResult{M,S,C}
    model::M
    solution::S
    config::C
    summaries::Dict{Symbol,Any}
    metadata::Dict{String,Any}
end
