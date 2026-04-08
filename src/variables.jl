"""
    VariableInfo(name, group, tags)

Metadata for a single state variable.

- `name`  — `Symbol` identifier used in output and layout queries
- `group` — `Symbol` grouping (e.g. `:mobile`, `:trap`); variables with the same
             group must be contiguous in the layout for `group_view` to work
- `tags`  — `Set{Symbol}` of capability tags, e.g. `:diffusion`, `:reaction`
"""
struct VariableInfo
    name::Symbol
    group::Symbol
    tags::Set{Symbol}
end

"""
    VariableLayout(variables)

Ordered collection of `VariableInfo` objects that defines how the flat state vector
maps onto `(nvariables, nx)`.  Groups must be contiguous.

Construct from a `Vector{VariableInfo}`; the group-range index is built automatically.
Query with `nvariables`, `variable_names`, `variables_in_group`, `variables_with_tag`.
"""
struct VariableLayout
    variables::Vector{VariableInfo}
    group_ranges::Dict{Symbol,UnitRange{Int}}
end

"""
    nvariables(layout) -> Int

Return the total number of variables in the layout.
"""
nvariables(layout::VariableLayout) = length(layout.variables)

"""
    variable_names(layout) -> Vector{Symbol}

Return an ordered vector of variable name symbols.
"""
variable_names(layout::VariableLayout) = [v.name for v in layout.variables]

"""
    has_group(layout, group) -> Bool

Return `true` if `group` is present in the layout.
"""
has_group(layout::VariableLayout, group::Symbol) = haskey(layout.group_ranges, group)

"""
    variables_in_group(layout, group) -> Vector{Int}

Return the variable indices (1-based) belonging to `group`.
Returns an empty vector if the group is not present.
"""
function variables_in_group(layout::VariableLayout, group::Symbol)
    has_group(layout, group) || return Int[]
    return collect(layout.group_ranges[group])
end

"""
    variables_with_tag(layout, tag) -> Vector{Int}

Return the variable indices (1-based) that carry `tag`.
"""
function variables_with_tag(layout::VariableLayout, tag::Symbol)
    idx = Int[]
    for (i, v) in enumerate(layout.variables)
        tag in v.tags && push!(idx, i)
    end
    return idx
end

function VariableLayout(variables::Vector{VariableInfo})
    group_ranges = Dict{Symbol,UnitRange{Int}}()

    if isempty(variables)
        return VariableLayout(variables, group_ranges)
    end

    start_idx = 1
    current_group = variables[1].group

    for i in 2:length(variables)
        if variables[i].group != current_group
            group_ranges[current_group] = start_idx:(i - 1)
            current_group = variables[i].group
            start_idx = i
        end
    end

    group_ranges[current_group] = start_idx:length(variables)

    return VariableLayout(variables, group_ranges)
end
