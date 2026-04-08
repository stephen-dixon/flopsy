using TOML

"""
    load_config(path) -> Dict

Load a TOML configuration file and return it as a `Dict`.
"""
function load_config(path::AbstractString)
    return TOML.parsefile(path)
end
