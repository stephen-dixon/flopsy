module Types

export SpeciesInfo, SimulationContext

struct SpeciesInfo
    names::Vector{String}
    nlevels::Int
end

struct SimulationContext
    species::SpeciesInfo
    nx::Int
    dx::Float64
end

end
