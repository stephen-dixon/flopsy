"""
    trapping_variable_layout(; mobile_name=:c, trap_name=:theta)

Construct a simple two-variable layout for a trapping-style model.

Variables:
- mobile variable: tagged with `:reaction` and `:diffusion`
- trap variable:   tagged with `:reaction`
"""
function trapping_variable_layout(; mobile_name::Symbol = :c, trap_name::Symbol = :theta)
    vars = [
        VariableInfo(mobile_name, :mobile, Set([:reaction, :diffusion])),
        VariableInfo(trap_name, :trap, Set([:reaction]))
    ]
    return VariableLayout(vars)
end

"""
    build_trapping_model(; ...)

Build a simple trapping-style reaction-diffusion model on top of the generic core.

This is intentionally a high-level convenience builder for a common model shape,
not a core concept of the framework.

Arguments:
- `mesh`: Mesh1D
- `k_trap`, `k_detrap`: reaction parameters
- `diffusion_coefficient`: diffusion coefficient applied to the mobile variable
- `mobile_name`, `trap_name`: symbolic variable names
- `aux`: additional auxiliary data
"""
function build_trapping_model(;
        mesh::Mesh1D,
        k_trap::Real,
        k_detrap::Real,
        diffusion_coefficient,
        mobile_name::Symbol = :c,
        trap_name::Symbol = :theta,
        aux = NamedTuple()
)
    layout = trapping_variable_layout(
        mobile_name = mobile_name,
        trap_name = trap_name
    )

    mobile_idx = 1
    trap_idx = 2

    reaction = SimpleTrappingReactionOperator(
        k_trap,
        k_detrap,
        mobile_idx,
        trap_idx
    )

    selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)

    coeffs = if diffusion_coefficient isa Real
        ConstantDiffusion([Float64(diffusion_coefficient), 0.0])
    elseif diffusion_coefficient isa AbstractDiffusionCoefficients
        diffusion_coefficient
    elseif diffusion_coefficient isa AbstractVector
        ConstantDiffusion(Float64.(diffusion_coefficient))
    else
        throw(ArgumentError("Unsupported diffusion_coefficient type $(typeof(diffusion_coefficient))"))
    end

    diffusion = LinearDiffusionOperator(coeffs, selector, nothing)

    return build_rd_model(
        layout = layout,
        mesh = mesh,
        reaction = reaction,
        diffusion = diffusion,
        aux = aux
    )
end
