using TOML

function load_config(path::AbstractString)
    return TOML.parsefile(path)
end
