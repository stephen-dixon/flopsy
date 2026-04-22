const _PLUGIN_PROVIDERS = Dict{Symbol, Function}()
const _CONFIG_DOMAINS = (:mesh, :backend, :ic, :bc, :output, :temperature, :problem)

struct ConfigValidationError <: Exception
    message::String
end

Base.showerror(io::IO, err::ConfigValidationError) = print(io, err.message)

"""
    register_syntax!(registry, spec)

Register a `SyntaxSpec` into the given registry.
"""
function register_syntax!(registry::SyntaxRegistry, spec::SyntaxSpec)
    key = (spec.domain, spec.type_name)
    registry.specs[key] = spec
    return registry
end

"""
    lookup_syntax(registry, domain, type_name)

Return the registered syntax entry for a `(domain, type)` pair.
"""
function lookup_syntax(registry::SyntaxRegistry, domain::Symbol, type_name::Symbol)
    key = (domain, type_name)
    haskey(registry.specs, key) ||
        throw(ConfigValidationError("Unknown syntax $(domain).$(type_name). Run `flopsy syntax list` to see registered syntax."))
    return registry.specs[key]
end

function syntax_specs(registry::SyntaxRegistry)
    return collect(values(registry.specs))
end

"""
    register_plugin_provider!(name, fn)

Register an in-process syntax provider, typically from a Julia package extension.
"""
function register_plugin_provider!(name::Symbol, fn::Function)
    _PLUGIN_PROVIDERS[name] = fn
    return name
end

function registered_plugin_providers()
    return collect(keys(_PLUGIN_PROVIDERS))
end

"""
    build_registry()

Construct a fresh syntax registry containing built-ins plus any loaded plugin
providers.
"""
function build_registry()
    registry = SyntaxRegistry()
    register_builtin_syntax!(registry)
    load_runtime_plugins!(registry)
    for name in sort!(collect(keys(_PLUGIN_PROVIDERS)); by = string)
        try
            _PLUGIN_PROVIDERS[name](registry)
        catch err
            throw(ConfigValidationError("Plugin provider `$(name)` failed during syntax registration: $(sprint(showerror, err))"))
        end
    end
    return registry
end

"""
    parse_input_deck(path)
    parse_input_deck(raw; path = nothing)

Parse a registry-driven TOML input deck into named blocks grouped by domain.
"""
function parse_input_deck(path::AbstractString)
    raw = TOML.parsefile(path)
    return parse_input_deck(raw; path = path)
end

function parse_input_deck(raw::AbstractDict; path = nothing)
    blocks = Dict{Symbol, Vector{ConfigBlock}}(domain => ConfigBlock[]
    for domain in _CONFIG_DOMAINS)

    for (key, value) in pairs(raw)
        domain = Symbol(key)
        domain in _CONFIG_DOMAINS ||
            throw(ConfigValidationError("Unknown top-level TOML domain '$key'. Supported domains: $(_CONFIG_DOMAINS)"))
        value isa AbstractDict ||
            throw(ConfigValidationError("Top-level domain '$key' must be a table of named blocks"))

        for (name, block) in pairs(value)
            block isa AbstractDict ||
                throw(ConfigValidationError("Block [$key.$name] must be a TOML table"))
            data = Dict{String, Any}(string(k) => v for (k, v) in pairs(block))
            haskey(data, "type") ||
                throw(ConfigValidationError("Block [$key.$name] is missing required field `type`"))
            data["type"] isa AbstractString ||
                throw(ConfigValidationError("Block [$key.$name] field `type` must be a string"))
            type_name = Symbol(data["type"])
            push!(blocks[domain], ConfigBlock(domain, Symbol(name), type_name, data))
        end
    end

    return InputDeck(path, Dict{String, Any}(string(k) => v for (k, v) in pairs(raw)), blocks)
end

"""
    build_context(deck; registry = build_registry())

Validate an `InputDeck` against the registry and build named runtime objects into
a `BuildContext`.
"""
function build_context(deck::InputDeck; registry::SyntaxRegistry = build_registry())
    ctx = BuildContext()
    for domain in _CONFIG_DOMAINS
        for block in get(deck.blocks, domain, ConfigBlock[])
            spec = lookup_syntax(registry, domain, block.type_name)
            data = _apply_parameter_defaults(spec, block.raw)
            _validate_block_schema(spec, data, block)
            _run_block_step(:validation, spec.validate, data, ctx, registry, block)
            obj = _run_block_step(:build, spec.build, data, ctx, registry, block)
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

function _validate_block_schema(spec::SyntaxSpec, data::Dict{String, Any}, block::ConfigBlock)
    _validate_unknown_fields(spec, data, block)
    for param in spec.parameters
        if param.required && !haskey(data, String(param.name))
            throw(ConfigValidationError(_format_block_field_message(
                block,
                param.name,
                "missing required field. Valid keys: $(_format_parameter_names(spec.parameters))"
            )))
        end
        haskey(data, String(param.name)) || continue
        _validate_parameter_value(param, data[String(param.name)], block, spec)
    end
end

function _validate_unknown_fields(spec::SyntaxSpec, data::Dict{String, Any}, block::ConfigBlock)
    allowed = Set(String(param.name) for param in spec.parameters)
    unknown = sort!(String[key for key in keys(data) if key ∉ allowed])
    isempty(unknown) && return nothing
    throw(ConfigValidationError(
        "Block [$(block.domain).$(block.name)] contains unknown key(s) $(join(repr.(unknown), ", ")). " *
        "Valid keys: $(_format_parameter_names(spec.parameters))",
    ))
end

function _validate_parameter_value(param::ParameterSpec, value, block::ConfigBlock, spec::SyntaxSpec)
    _value_matches_kind(param.kind, value, param.element_kind) ||
        throw(ConfigValidationError(_format_block_field_message(
            block,
            param.name,
            "expected $(_describe_kind(param.kind, param.element_kind)), got $(typeof(value))"
        )))

    if param.allowed_values !== nothing
        normalized = _normalize_allowed_value(value)
        allowed = [_normalize_allowed_value(v) for v in param.allowed_values]
        normalized in allowed || throw(ConfigValidationError(_format_block_field_message(
            block,
            param.name,
            "invalid value $(repr(value)). Allowed values: $(join(repr.(param.allowed_values), ", "))"
        )))
    end

    return nothing
end

function _run_block_step(
        stage::Symbol, fn::Function, data, ctx, registry, block::ConfigBlock)
    try
        return fn(data, ctx, registry, block)
    catch err
        err isa ConfigValidationError && rethrow()
        throw(ConfigValidationError(
            "Block [$(block.domain).$(block.name)] failed during $(stage): $(sprint(showerror, err))",
        ))
    end
end

function _format_parameter_names(params::Vector{ParameterSpec})
    return join(sort!(string.(getfield.(params, :name))), ", ")
end

function _format_block_field_message(block::ConfigBlock, field::Symbol, message::AbstractString)
    return "Block [$(block.domain).$(block.name)] field `$(field)` $message"
end

function _value_matches_kind(kind::Symbol, value, element_kind::Symbol)
    if kind == :any
        return true
    elseif kind == :string
        return value isa AbstractString
    elseif kind == :integer
        return value isa Integer
    elseif kind == :real
        return value isa Real
    elseif kind == :bool
        return value isa Bool
    elseif kind == :vector
        value isa AbstractVector || return false
        return all(item -> _value_matches_kind(element_kind, item, :any), value)
    end
    return true
end

function _describe_kind(kind::Symbol, element_kind::Symbol)
    if kind == :vector
        return "vector of $(_describe_kind(element_kind, :any)) values"
    end
    return String(kind)
end

_normalize_allowed_value(value) = value isa AbstractString ? String(value) : value

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
    elseif domain == :temperature
        ctx.temperatures
    elseif domain == :problem
        ctx.problems
    else
        error("Unsupported build-context domain $domain")
    end
    haskey(storage, name) &&
        throw(ConfigValidationError("Duplicate object name '$name' in domain '$domain'"))
    storage[name] = obj
    return obj
end

"""
    syntax_list(registry = build_registry())

Return a sorted list of registered `(domain, type, plugin)` tuples.
"""
function syntax_list(registry::SyntaxRegistry = build_registry())
    rows = [(spec.domain, spec.type_name, spec.plugin) for spec in syntax_specs(registry)]
    sort!(rows; by = x -> (string(x[1]), string(x[2])))
    return rows
end

"""
    syntax_show(domain, type_name; registry = build_registry())

Return detailed information for one registered syntax entry.
"""
function syntax_show(domain::Symbol, type_name::Symbol; registry::SyntaxRegistry = build_registry())
    spec = lookup_syntax(registry, domain, type_name)
    return (
        domain = spec.domain,
        type_name = spec.type_name,
        plugin = spec.plugin,
        help = spec.help,
        parameters = spec.parameters
    )
end
