function build_hotgates_trapping_model(;
    mesh::Mesh1D,
    model,
    adaptor::HotgatesTrappingAdaptor,
    temperature::AbstractTemperatureProvider,
    diffusion_coefficients::AbstractVector,
)
    layout = build_variable_layout(adaptor)

    reaction = HotgatesReactionOperator(model, adaptor, temperature)

    selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)

    diffusion = LinearDiffusionOperator(diffusion_coefficients, selector)

    return build_rd_model(
        layout = layout,
        mesh = mesh,
        reaction = reaction,
        diffusion = diffusion,
    )
end
