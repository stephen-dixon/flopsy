const _PLUGIN_PROVIDERS = Dict{Symbol, Function}()
const _CONFIG_DOMAINS = (:mesh, :backend, :ic, :bc, :output, :problem)

function register_syntax!(registry::SyntaxRegistry, spec::SyntaxSpec)
    key = (spec.domain, spec.type_name)
    registry.specs[key] = spec
    return registry
end

function lookup_syntax(registry::SyntaxRegistry, domain::Symbol, type_name::Symbol)
    key = (domain, type_name)
    haskey(registry.specs, key) ||
        throw(ArgumentError("Unknown syntax $(domain).$(type_name). Run `flopsy syntax list` to see registered syntax."))
    return registry.specs[key]
end

function syntax_specs(registry::SyntaxRegistry)
    return collect(values(registry.specs))
end

function register_plugin_provider!(name::Symbol, fn::Function)
    _PLUGIN_PROVIDERS[name] = fn
    return name
end

function registered_plugin_providers()
    return collect(keys(_PLUGIN_PROVIDERS))
end

function build_registry()
    registry = SyntaxRegistry()
    register_builtin_syntax!(registry)
    load_runtime_plugins!(registry)
    for name in sort!(collect(keys(_PLUGIN_PROVIDERS)); by = string)
        _PLUGIN_PROVIDERS[name](registry)
    end
    return registry
end

function parse_input_deck(path::AbstractString)
    raw = TOML.parsefile(path)
    return parse_input_deck(raw; path = path)
end

function parse_input_deck(raw::AbstractDict; path = nothing)
    blocks = Dict{Symbol, Vector{ConfigBlock}}(domain => ConfigBlock[] for domain in _CONFIG_DOMAINS)

    for (key, value) in pairs(raw)
        domain = Symbol(key)
        domain in _CONFIG_DOMAINS ||
            throw(ArgumentError("Unknown top-level TOML domain '$key'. Supported domains: $(_CONFIG_DOMAINS)"))
        value isa AbstractDict ||
            throw(ArgumentError("Top-level domain '$key' must be a table of named blocks"))

        for (name, block) in pairs(value)
            block isa AbstractDict ||
                throw(ArgumentError("Block [$key.$name] must be a TOML table"))
            data = Dict{String, Any}(string(k) => v for (k, v) in pairs(block))
            haskey(data, "type") ||
                throw(ArgumentError("Block [$key.$name] is missing required field `type`"))
            type_name = Symbol(data["type"])
            push!(blocks[domain], ConfigBlock(domain, Symbol(name), type_name, data))
        end
    end

    return InputDeck(path, Dict{String, Any}(string(k) => v for (k, v) in pairs(raw)), blocks)
end

function build_context(deck::InputDeck; registry::SyntaxRegistry = build_registry())
    ctx = BuildContext()
    for domain in _CONFIG_DOMAINS
        for block in get(deck.blocks, domain, ConfigBlock[])
            spec = lookup_syntax(registry, domain, block.type_name)
            data = _apply_parameter_defaults(spec, block.raw)
            _validate_required_fields(spec, data, block)
            spec.validate(data, ctx, registry, block)
            obj = spec.build(data, ctx, registry, block)
            _store_built_object!(ctx, domain, block.name, obj)
        end
    end
    ctx.artifacts[:deck] = deck
    ctx.artifacts[:registry] = registry
    return ctx
end

function _apply_parameter_defaults(spec::SyntaxSpec, raw::Dict{String, Any})
    data = copy(raw)
    for param in spec.parameters
        if !haskey(data, String(param.name)) && param.default !== nothing
            data[String(param.name)] = param.default
        end
    end
    return data
end

function _validate_required_fields(spec::SyntaxSpec, data::Dict{String, Any}, block::ConfigBlock)
    for param in spec.parameters
        if param.required && !haskey(data, String(param.name))
            throw(ArgumentError("Block [$(block.domain).$(block.name)] missing required field `$(param.name)`"))
        end
    end
end

function _store_built_object!(ctx::BuildContext, domain::Symbol, name::Symbol, obj)
    storage = if domain == :mesh
        ctx.meshes
    elseif domain == :backend
        ctx.backends
    elseif domain == :ic
        ctx.ics
    elseif domain == :bc
        ctx.bcs
    elseif domain == :output
        ctx.outputs
    elseif domain == :problem
        ctx.problems
    else
        error("Unsupported build-context domain $domain")
    end
    haskey(storage, name) &&
        throw(ArgumentError("Duplicate object name '$name' in domain '$domain'"))
    storage[name] = obj
    return obj
end

function syntax_list(registry::SyntaxRegistry = build_registry())
    rows = [(spec.domain, spec.type_name, spec.plugin) for spec in syntax_specs(registry)]
    sort!(rows; by = x -> (string(x[1]), string(x[2])))
    return rows
end

function syntax_show(domain::Symbol, type_name::Symbol; registry::SyntaxRegistry = build_registry())
    spec = lookup_syntax(registry, domain, type_name)
    return (
        domain = spec.domain,
        type_name = spec.type_name,
        plugin = spec.plugin,
        help = spec.help,
        parameters = spec.parameters,
    )
end
