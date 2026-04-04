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

struct SystemContext{T, M, A}
    nx::Int
    mesh::M
    aux::A
    scratch::Dict{Symbol, Any}
end

struct SystemModel{L, O, C} <: AbstractSystemModel
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
    kwargs::Dict{Symbol, Any} = Dict{Symbol, Any}()
end
