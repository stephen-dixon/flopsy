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
struct InitialConditionConfig{T}
    kind::Symbol
    amplitude::Union{Nothing, T}
    value::Union{Nothing, T}
    mobile_amplitude::Union{Nothing, T}
    mobile_value::Union{Nothing, T}
    trap_occupancy::Union{Nothing, T}
end

function InitialConditionConfig(
    kind::Symbol,
    amplitude::Union{Nothing, Real},
    value::Union{Nothing, Real},
    mobile_amplitude::Union{Nothing, Real},
    mobile_value::Union{Nothing, Real},
    trap_occupancy::Union{Nothing, Real})
    vals = (
        amplitude,
        value,
        mobile_amplitude,
        mobile_value,
        trap_occupancy,
    )
    T = promote_type(
        map(v -> v === nothing ? Float64 : typeof(v), vals)...,
    )
    return InitialConditionConfig{T}(
        kind,
        amplitude,
        value,
        mobile_amplitude,
        mobile_value,
        trap_occupancy,
    )
end

struct ProblemConfig{M, S, BC, IC, P}
    problem_type::Symbol
    mesh::M
    solver::S
    boundary_conditions::BC
    initial_conditions::IC
    parameters::P
end
