function build_rd_model(;
    layout::VariableLayout,
    mesh::Mesh1D,
    reaction::Union{AbstractOperator,Nothing}=nothing,
    diffusion::Union{AbstractOperator,Nothing}=nothing,
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
        reaction = reaction,
        diffusion = diffusion,
        constraints = constraints,
    )

    return SystemModel(layout, operators, ctx)
end
