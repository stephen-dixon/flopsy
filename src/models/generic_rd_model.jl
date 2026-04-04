function build_rd_model(;
    layout::VariableLayout,
    mesh::Mesh1D,
    reaction::Union{AbstractOperator, Nothing}=nothing,
    diffusion::Union{AbstractOperator, Nothing}=nothing,
    constraints::Union{AbstractOperator, Nothing}=nothing,
    aux::Dict{Symbol, Any}=Dict{Symbol, Any}(),
)
    nx = length(mesh.x)

    aux2 = copy(aux)
    aux2[:layout] = layout

    ctx = SystemContext(
        nx,
        mesh,
        aux2,
        Dict{Symbol, Any}(),
    )

    operators = (
        reaction = reaction,
        diffusion = diffusion,
        constraints = constraints,
    )

    return SystemModel(layout, operators, ctx)
end
