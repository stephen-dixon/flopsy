struct ParameterSpec
    name::Symbol
    required::Bool
    default
    doc::String
end

struct SyntaxSpec{V, B}
    domain::Symbol
    type_name::Symbol
    parameters::Vector{ParameterSpec}
    help::String
    validate::V
    build::B
    plugin::Symbol
end

struct SyntaxRegistry
    specs::Dict{Tuple{Symbol, Symbol}, SyntaxSpec}
end

SyntaxRegistry() = SyntaxRegistry(Dict{Tuple{Symbol, Symbol}, SyntaxSpec}())

struct ConfigBlock
    domain::Symbol
    name::Symbol
    type_name::Symbol
    raw::Dict{String, Any}
end

struct InputDeck
    path::Union{Nothing, String}
    raw::Dict{String, Any}
    blocks::Dict{Symbol, Vector{ConfigBlock}}
end

mutable struct BuildContext
    meshes::Dict{Symbol, Any}
    backends::Dict{Symbol, Any}
    ics::Dict{Symbol, Any}
    bcs::Dict{Symbol, Any}
    outputs::Dict{Symbol, Any}
    problems::Dict{Symbol, Any}
    artifacts::Dict{Symbol, Any}
end

BuildContext() = BuildContext(
    Dict{Symbol, Any}(),
    Dict{Symbol, Any}(),
    Dict{Symbol, Any}(),
    Dict{Symbol, Any}(),
    Dict{Symbol, Any}(),
    Dict{Symbol, Any}(),
    Dict{Symbol, Any}(),
)

struct SpeciesInfo
    name::Symbol
    role::Symbol
    transport::Symbol
    index::Int
    doc::String
    boundary_target::Bool
end

struct BackendDefinition{S, F, M}
    name::Symbol
    type_name::Symbol
    species::S
    build_model::F
    metadata::M
end

struct ICDefinition{A, F}
    name::Symbol
    type_name::Symbol
    backend_ref::Union{Nothing, Symbol}
    affects::A
    apply!::F
end

struct BCDefinition
    name::Symbol
    type_name::Symbol
    species::Symbol
    boundary::Symbol
    value::Float64
    method::Symbol
end

struct OutputDefinition
    name::Symbol
    type_name::Symbol
    file::String
    xdmf::Bool
end

struct ProblemDefinition{S}
    name::Symbol
    type_name::Symbol
    mesh::Symbol
    backend::Symbol
    ics::Vector{Symbol}
    bcs::Vector{Symbol}
    outputs::Vector{Symbol}
    tspan::Tuple{Float64, Float64}
    solver::S
end

struct ConfiguredSimulation{P, O, C}
    problem::P
    outputs::O
    context::C
end
