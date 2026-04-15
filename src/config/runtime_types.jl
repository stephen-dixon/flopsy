"""
    ParameterSpec(name, required, default, doc; kind = :any, allowed_values = nothing, element_kind = :any)

Declarative schema for one user-facing TOML field in a registered syntax block.
"""
struct ParameterSpec
    name::Symbol
    required::Bool
    default
    doc::String
    kind::Symbol
    allowed_values::Union{Nothing, Vector{Any}}
    element_kind::Symbol
end

function ParameterSpec(
    name::Symbol,
    required::Bool,
    default,
    doc::AbstractString;
    kind::Symbol = :any,
    allowed_values = nothing,
    element_kind::Symbol = :any,
)
    allowed = allowed_values === nothing ? nothing : collect(allowed_values)
    return ParameterSpec(name, required, default, String(doc), kind, allowed, element_kind)
end

"""
    SyntaxSpec

Registry entry describing a supported `(domain, type)` input-deck block.
"""
struct SyntaxSpec{V, B}
    domain::Symbol
    type_name::Symbol
    parameters::Vector{ParameterSpec}
    help::String
    validate::V
    build::B
    plugin::Symbol
end

"""
    SyntaxRegistry()

In-memory registry of built-in and plugin-provided input-deck syntax.
"""
struct SyntaxRegistry
    specs::Dict{Tuple{Symbol, Symbol}, SyntaxSpec}
end

SyntaxRegistry() = SyntaxRegistry(Dict{Tuple{Symbol, Symbol}, SyntaxSpec}())

"""
    ConfigBlock

Parsed named TOML block before validation/build.
"""
struct ConfigBlock
    domain::Symbol
    name::Symbol
    type_name::Symbol
    raw::Dict{String, Any}
end

"""
    InputDeck

Parsed registry-driven TOML input deck grouped by top-level domain.
"""
struct InputDeck
    path::Union{Nothing, String}
    raw::Dict{String, Any}
    blocks::Dict{Symbol, Vector{ConfigBlock}}
end

"""
    BuildContext()

Validated named object graph produced from an `InputDeck`.
"""
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

"""
    BackendDefinition

Named backend definition produced from a registered backend syntax block.
"""
struct BackendDefinition{S, F, M}
    name::Symbol
    type_name::Symbol
    species::S
    build_model::F
    metadata::M
end

"""
    ICDefinition

Named initial-condition definition produced from a registered IC syntax block.
"""
struct ICDefinition{A, F}
    name::Symbol
    type_name::Symbol
    backend_ref::Union{Nothing, Symbol}
    affects::A
    apply!::F
end

"""
    BCDefinition

Named boundary-condition definition produced from a registered BC syntax block.
"""
struct BCDefinition
    name::Symbol
    type_name::Symbol
    species::Symbol
    boundary::Symbol
    value::Float64
    method::Symbol
end

"""
    OutputDefinition

Named output definition produced from a registered output syntax block.
"""
struct OutputDefinition
    name::Symbol
    type_name::Symbol
    file::String
    xdmf::Bool
end

"""
    ProblemDefinition

Named problem assembly block referencing previously built objects.
"""
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

"""
    ConfiguredSimulation

Assembled simulation plan consisting of a `SimulationProblem`, requested outputs,
and the originating build context.
"""
struct ConfiguredSimulation{P, O, C}
    problem::P
    outputs::O
    context::C
end
