module RHS

using ..Reaction
using ..Diffusion
using ..State

export rhs!

function rhs!(du, u, p, t)
    ctx = p

    fill!(du, 0.0)

    # Reaction contribution
    Reaction.reaction_rates!(du, u, ctx, t)

    # Diffusion contribution
    Diffusion.diffusion_term!(du, u, ctx)
end

end
