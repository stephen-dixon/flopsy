"""
    state_view(u, layout, nx) -> Matrix

Reshape the flat state vector `u` (length `nvariables * nx`) into a
`(nvariables, nx)` matrix without copying.  Column `ix` is the state at node `ix`.
"""
function state_view(u::AbstractVector, layout::VariableLayout, nx::Integer)
    reshape(u, nvariables(layout), nx)
end

"""
    node_view(U, ix) -> AbstractVector

Return a view of the state matrix `U` at node index `ix` (column `ix`).
"""
node_view(U::AbstractMatrix, ix::Integer) = @view U[:, ix]

"""
    variable_view(U, ivar) -> AbstractVector

Return a view of the state matrix `U` for variable index `ivar` across all nodes
(row `ivar`).
"""
variable_view(U::AbstractMatrix, ivar::Integer) = @view U[ivar, :]

"""
    group_view(U, layout, group) -> AbstractMatrix

Return a view of the state matrix `U` containing only the rows belonging to `group`.
"""
function group_view(U::AbstractMatrix, layout::VariableLayout, group::Symbol)
    r = layout.group_ranges[group]
    return @view U[r, :]
end
