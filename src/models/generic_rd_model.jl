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
    aux::Dict{Symbol,Any}=Dict{Symbol,Any}(),
)
    nx = length(mesh.x)

    ctx = SystemContext(
        layout,
        nx,
        mesh,
        copy(aux),
        Dict{Symbol,Any}(),
    )

    operators = (
        reaction    = reaction,
        diffusion   = diffusion,
        boundary    = boundary,
        constraints = constraints,
    )

    return SystemModel(layout, operators, ctx)
end
