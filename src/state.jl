module State

export create_initial_state, reshape_state

"""
State layout:
flattened vector for solver
but logically: (species, level, node)
"""

function create_initial_state(nspecies, nlevels, nx)
    return zeros(nspecies * nlevels * nx)
end

function reshape_state(u, nspecies, nlevels, nx)
    return reshape(u, nspecies, nlevels, nx)
end

end
