struct VariableInfo
    name::Symbol
    group::Symbol
    tags::Set{Symbol}
end

struct VariableLayout
    variables::Vector{VariableInfo}
    group_ranges::Dict{Symbol,UnitRange{Int}}
end

nvariables(layout::VariableLayout) = length(layout.variables)

variable_names(layout::VariableLayout) = [v.name for v in layout.variables]

has_group(layout::VariableLayout, group::Symbol) = haskey(layout.group_ranges, group)

function variables_in_group(layout::VariableLayout, group::Symbol)
    has_group(layout, group) || return Int[]
    return collect(layout.group_ranges[group])
end

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
