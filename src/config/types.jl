"""
    MeshConfig

Legacy typed mesh specification for the deprecated config-driven problem layer.
"""
struct MeshConfig{T}
    kind::Symbol
    nx::Int
    xmin::T
    xmax::T
end

"""
    BoundaryConditionConfig

Legacy typed boundary-condition specification used by the deprecated config
layer.
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

Legacy typed solver settings parsed from TOML before they are lowered into the
core `SolverConfig`.
"""
struct InputSolverConfig{A, T, S}
    formulation::Symbol
    algorithm::A
    dt::Union{Nothing, T}
    abstol::T
    reltol::T
    saveat::S
    split_method::Symbol
    show_progress::Bool
end

# Backward-compatible constructor without split_method (defaults to Strang splitting)
function InputSolverConfig(formulation, algorithm, dt, abstol, reltol, saveat)
    return InputSolverConfig(
        formulation, algorithm, dt, abstol, reltol, saveat, :strang, true)
end

function InputSolverConfig(formulation, algorithm, dt, abstol, reltol, saveat, split_method)
    return InputSolverConfig(
        formulation, algorithm, dt, abstol, reltol, saveat, split_method, true)
end

"""
    InitialConditionConfig

Legacy typed initial-condition specification used by the deprecated factory
built problem layer.
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
        trap_occupancy
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
        trap_occupancy
    )
end

"""
    ProblemConfig

Legacy typed top-level problem configuration used by the deprecated factory
built problem layer.
"""
struct ProblemConfig{M, S, BC, IC, P}
    problem_type::Symbol
    mesh::M
    solver::S
    boundary_conditions::BC
    initial_conditions::IC
    parameters::P
end
