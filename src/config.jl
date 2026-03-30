module Config

using TOML

export SimulationConfig, load_config

struct SimulationConfig
    tspan::Tuple{Float64, Float64}
    nx::Int
    length::Float64
    solver::String
    abstol::Float64
    reltol::Float64
    output_stride::Int
    reaction_lib::String
end

function load_config(path::String)
    data = TOML.parsefile(path)

    return SimulationConfig(
        (data["time"]["t0"], data["time"]["tend"]),
        data["mesh"]["nx"],
        data["mesh"]["length"],
        data["solver"]["method"],
        data["solver"]["abstol"],
        data["solver"]["reltol"],
        data["output"]["stride"],
        data["reaction"]["library"]
    )
end

end
