"""
    MeshConfig

Typed mesh specification for the config-driven problem layer.
"""
struct MeshConfig{T}
    kind::Symbol
    nx::Int
    xmin::T
    xmax::T
end

"""
    BoundaryConditionConfig

Typed boundary-condition specification used by the config layer.
"""
struct BoundaryConditionConfig{T}
    variable::Symbol
    side::Symbol
    kind::Symbol
    value::T
    method::Symbol
end

"""
    InputSolverConfig

Typed solver settings parsed from TOML before they are lowered into the core
`SolverConfig`.
"""
struct InputSolverConfig{A, T, S}
    formulation::Symbol
    algorithm::A
    dt::Union{Nothing, T}
    abstol::T
    reltol::T
    saveat::S
end

"""
    ProblemConfig

Typed top-level problem configuration used by the factory-built problem layer.
"""
struct ProblemConfig{M, S, BC, P}
    problem_type::Symbol
    mesh::M
    solver::S
    boundary_conditions::BC
    parameters::P
end
