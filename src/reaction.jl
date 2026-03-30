module Reaction

export reaction_rates!

"""
Placeholder: replace with ccall into external library
"""

function reaction_rates!(du, u, ctx, t)
    # Example: no-op
    @inbounds for i in eachindex(u)
        du[i] += 0.0
    end
end

end
