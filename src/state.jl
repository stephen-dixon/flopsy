function state_view(u::AbstractVector, layout::VariableLayout, nx::Integer)
    reshape(u, nvariables(layout), nx)
end

node_view(U::AbstractMatrix, ix::Integer) = @view U[:, ix]
variable_view(U::AbstractMatrix, ivar::Integer) = @view U[ivar, :]

function group_view(U::AbstractMatrix, layout::VariableLayout, group::Symbol)
    r = layout.group_ranges[group]
    return @view U[r, :]
end
