"""
    build_rd_model(; layout, mesh, reaction=nothing, diffusion=nothing,
                     boundary=nothing, constraints=nothing, aux=Dict())

Construct a `SystemModel` from a variable layout, mesh, and operator slots.

Each slot may be `nothing` (omitted) or any `AbstractOperator`.  The boundary slot
is intended for a `DirichletBoundaryOperator`; the constraints slot for a
`ConstraintOperator` (DAE path).
"""
function build_rd_model(;
    layout::VariableLayout,
    mesh::Mesh1D,
    reaction::Union{AbstractOperator,Nothing}=nothing,
    diffusion::Union{AbstractOperator,Nothing}=nothing,
    boundary::Union{AbstractOperator,Nothing}=nothing,
    constraints::Union{AbstractOperator,Nothing}=nothing,
    aux=NamedTuple(),
    scratch=nothing,
)
    nx = length(mesh.x)
    nstate = nvariables(layout) * nx
    default_scratch = (
        rhs_tmp = zeros(Float64, nstate),
        implicit_rhs_tmp = zeros(Float64, nstate),
    )
    scratch_data = if scratch === nothing
        default_scratch
    elseif scratch isa NamedTuple
        merge(default_scratch, scratch)
    else
        scratch
    end

    ctx = SystemContext(
        layout,
        nx,
        mesh,
        aux isa AbstractDict ? copy(aux) : aux,
        scratch_data,
    )

    operators = (
        reaction    = reaction,
        diffusion   = diffusion,
        boundary    = boundary,
        constraints = constraints,
    )

    return SystemModel(layout, operators, ctx)
end
