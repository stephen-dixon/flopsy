struct MeshConfig{T}
    kind::Symbol
    nx::Int
    xmin::T
    xmax::T
end

struct BoundaryConditionConfig{T}
    variable::Symbol
    side::Symbol
    kind::Symbol
    value::T
    method::Symbol
end

struct InputSolverConfig{A,T,S}
    formulation::Symbol
    algorithm::A
    dt::Union{Nothing,T}
    abstol::T
    reltol::T
    saveat::S
end

struct ProblemConfig{M,S,BC,P}
    problem_type::Symbol
    mesh::M
    solver::S
    boundary_conditions::BC
    parameters::P
end
