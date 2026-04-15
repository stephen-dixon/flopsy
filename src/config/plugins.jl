using Pkg

function plugin_home()
    return joinpath(homedir(), ".flopsy")
end

function plugin_environment_path()
    return joinpath(plugin_home(), "plugins")
end

function plugin_state_file()
    return joinpath(plugin_environment_path(), "plugins.toml")
end

function ensure_plugin_environment()
    env = plugin_environment_path()
    mkpath(env)
    project = joinpath(env, "Project.toml")
    if !isfile(project)
        write(project, "name = \"FlopsyRuntimePlugins\"\nuuid = \"cb90c8ef-09af-4e3a-b1b5-d21c83d31990\"\nversion = \"0.1.0\"\n")
    end
    state = plugin_state_file()
    if !isfile(state)
        write(state, "[plugins]\n")
    end
    return env
end

function _with_plugin_env(f::Function)
    env = ensure_plugin_environment()
    prev = Base.active_project()
    try
        Pkg.activate(env; io = devnull)
        return f(env)
    finally
        prev === nothing ? nothing : Pkg.activate(dirname(prev); io = devnull)
    end
end

function _read_plugin_state()
    ensure_plugin_environment()
    return TOML.parsefile(plugin_state_file())
end

function _write_plugin_state(state::Dict{String, Any})
    open(plugin_state_file(), "w") do io
        TOML.print(io, state)
    end
end

"""
    plugin_list()

List runtime plugins known to the managed Flopsy plugin environment.
"""
function plugin_list()
    state = _read_plugin_state()
    plugins = get(state, "plugins", Dict{String, Any}())
    out = NamedTuple[]
    for name in sort!(collect(keys(plugins)))
        info = plugins[name]
        push!(out, (
            name = name,
            source = get(info, "source", ""),
            registry = get(info, "registry", ""),
            load_error = get(info, "load_error", ""),
        ))
    end
    return out
end

"""
    plugin_register!(name; registry = nothing, url = nothing, path = nothing)

Register and install a runtime plugin into the managed Flopsy plugin environment.
"""
function plugin_register!(name::AbstractString; registry = nothing, url = nothing, path = nothing)
    path !== nothing && !isdir(path) &&
        throw(ArgumentError("Plugin path does not exist or is not a directory: $(path)"))
    _with_plugin_env() do _
        registry !== nothing && Pkg.Registry.add(Pkg.RegistrySpec(url = registry); io = devnull)
        if path !== nothing
            Pkg.develop(Pkg.PackageSpec(path = path); io = devnull)
        elseif url !== nothing
            Pkg.add(Pkg.PackageSpec(url = url); io = devnull)
        else
            Pkg.add(Pkg.PackageSpec(name = String(name)); io = devnull)
        end
    end

    state = _read_plugin_state()
    plugins = get!(state, "plugins") do
        Dict{String, Any}()
    end
    plugins[String(name)] = Dict(
        "source" => path !== nothing ? path : url !== nothing ? url : String(name),
        "registry" => registry === nothing ? "" : String(registry),
        "load_error" => "",
    )
    _write_plugin_state(state)
    return name
end

"""
    plugin_remove!(name)

Remove a runtime plugin from the managed Flopsy plugin environment.
"""
function plugin_remove!(name::AbstractString)
    _with_plugin_env() do _
        Pkg.rm(Pkg.PackageSpec(name = String(name)); io = devnull)
    end
    state = _read_plugin_state()
    plugins = get(state, "plugins", Dict{String, Any}())
    delete!(plugins, String(name))
    state["plugins"] = plugins
    _write_plugin_state(state)
    return name
end

function load_runtime_plugins!(registry::SyntaxRegistry)
    state = _read_plugin_state()
    plugins = get(state, "plugins", Dict{String, Any}())
    isempty(plugins) && return registry

    env = ensure_plugin_environment()
    env in LOAD_PATH || pushfirst!(LOAD_PATH, env)

    dirty = false
    for name in keys(plugins)
        info = plugins[name]
        try
            mod = _import_runtime_plugin(name)
            if isdefined(mod, :register_flopsy_plugin!)
                getfield(mod, :register_flopsy_plugin!)(registry)
            else
                throw(ConfigValidationError("Plugin `$name` does not define `register_flopsy_plugin!(registry)`"))
            end
            info["load_error"] = ""
        catch err
            info["load_error"] = "plugin `$name` failed to load: " * sprint(showerror, err)
            dirty = true
        end
    end
    dirty && _write_plugin_state(state)
    return registry
end

function _import_runtime_plugin(name::AbstractString)
    sym = Symbol(name)
    Core.eval(@__MODULE__, Expr(:import, sym))
    return getfield(@__MODULE__, sym)
end
