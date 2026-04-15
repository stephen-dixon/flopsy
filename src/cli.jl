const _CLI_MAIN_HELP = """
flopsy: registry-driven reaction-diffusion input-deck runner

Usage:
  flopsy run <input.toml>
  flopsy validate <input.toml>
  flopsy syntax list
  flopsy syntax show <domain> <type>
  flopsy plugin list
  flopsy plugin register <name> [--registry <url>] [--url <pkg-url> | --path <local-path>]
  flopsy plugin remove <name>
  flopsy xdmf <fields.h5> [--output out.xdmf]
  flopsy --help
  flopsy --version

Main workflows:
  run       Validate an input deck and execute the selected [problem.<name>] block.
  validate  Parse and assemble the input deck without running the solver.
  syntax    Inspect built-in and plugin-provided registered syntax.
  plugin    Manage runtime plugins loaded into the syntax registry.
  xdmf      Generate an XDMF companion for an HDF5 field output file.

Use `flopsy <command> --help` for command-specific help.
"""

const _CLI_COMMAND_HELP = Dict(
    "run" => """
Usage:
  flopsy run <input.toml>

Run a registry-driven input deck from a TOML file.
""",
    "validate" => """
Usage:
  flopsy validate <input.toml>

Validate a registry-driven input deck without running the solver.
""",
    "syntax" => """
Usage:
  flopsy syntax list
  flopsy syntax show <domain> <type>

Inspect the registered input-deck syntax. Output reflects built-ins plus any
successfully loaded runtime plugins.
""",
    "plugin" => """
Usage:
  flopsy plugin list
  flopsy plugin register <name> [--registry <url>] [--url <pkg-url> | --path <local-path>]
  flopsy plugin remove <name>

Plugins are installed into a managed runtime Julia environment under the user
config directory and can register additional syntax into the live registry.
""",
    "xdmf" => """
Usage:
  flopsy xdmf <fields.h5> [--output out.xdmf]

Generate an XDMF companion file for a Flopsy HDF5 field output.
""",
)

"""
    cli_main(args = copy(ARGS))

Run the Flopsy command-line interface and return a process-style exit code.
"""
function cli_main(args::Vector{String} = copy(ARGS))
    verbose = false
    filtered = String[]
    for arg in args
        if arg == "--verbose"
            verbose = true
        else
            push!(filtered, arg)
        end
    end

    try
        return _cli_dispatch(filtered)
    catch err
        if verbose
            showerror(stderr, err, catch_backtrace())
            println(stderr)
        else
            println(stderr, "error: ", sprint(showerror, err))
            println(stderr, "hint: rerun with --verbose for a stack trace")
        end
        return 1
    end
end

function _cli_dispatch(args::Vector{String})
    isempty(args) && return _cli_usage(0)
    cmd = popfirst!(args)

    cmd in ("--help", "-h", "help") && return _cli_usage(0)
    cmd == "--version" && return _cli_version()

    if cmd == "run"
        return _cli_run(args)
    elseif cmd == "validate"
        return _cli_validate(args)
    elseif cmd == "syntax"
        return _cli_syntax(args)
    elseif cmd == "plugin"
        return _cli_plugin(args)
    elseif cmd == "xdmf"
        return _cli_xdmf(args)
    end

    throw(ArgumentError("Unknown command `$cmd`. Use `flopsy --help` to see available commands."))
end

function _cli_usage(code::Int)
    println(_CLI_MAIN_HELP)
    return code
end

function _cli_command_usage(cmd::AbstractString)
    println(get(_CLI_COMMAND_HELP, String(cmd), _CLI_MAIN_HELP))
    return 0
end

function _cli_version()
    version = something(Base.pkgversion(@__MODULE__), v"0.0.0")
    println("flopsy ", version)
    return 0
end

function _cli_run(args::Vector{String})
    _wants_help(args) && return _cli_command_usage("run")
    length(args) == 1 || throw(ArgumentError("`flopsy run` expects exactly one input deck path."))
    run_input_deck(args[1])
    return 0
end

function _cli_validate(args::Vector{String})
    _wants_help(args) && return _cli_command_usage("validate")
    length(args) == 1 || throw(ArgumentError("`flopsy validate` expects exactly one input deck path."))
    validate_input_deck(args[1])
    println("valid")
    return 0
end

function _cli_syntax(args::Vector{String})
    isempty(args) && return _cli_command_usage("syntax")
    _wants_help(args) && return _cli_command_usage("syntax")

    sub = popfirst!(args)
    if sub == "list"
        isempty(args) || throw(ArgumentError("`flopsy syntax list` does not accept extra arguments."))
        for row in syntax_list()
            println("$(row[1])\t$(row[2])\tplugin=$(row[3])")
        end
        return 0
    elseif sub == "show"
        length(args) == 2 || throw(ArgumentError("`flopsy syntax show` expects <domain> <type>."))
        info = syntax_show(Symbol(args[1]), Symbol(args[2]))
        println("domain: ", info.domain)
        println("type: ", info.type_name)
        println("plugin: ", info.plugin)
        println("help: ", info.help)
        println("parameters:")
        for param in info.parameters
            print("  - ", param.name)
            print(" required=", param.required)
            print(" default=", repr(param.default))
            print(" kind=", param.kind)
            if param.kind == :vector
                print("[", param.element_kind, "]")
            end
            param.allowed_values === nothing || print(" allowed=", repr(param.allowed_values))
            println(" :: ", param.doc)
        end
        return 0
    end

    throw(ArgumentError("Unknown syntax subcommand `$sub`. Use `flopsy syntax --help`."))
end

function _cli_plugin(args::Vector{String})
    isempty(args) && return _cli_command_usage("plugin")
    _wants_help(args) && return _cli_command_usage("plugin")

    sub = popfirst!(args)
    if sub == "list"
        isempty(args) || throw(ArgumentError("`flopsy plugin list` does not accept extra arguments."))
        for info in plugin_list()
            line = info.name * "\tsource=" * info.source
            !isempty(info.registry) && (line *= "\tregistry=" * info.registry)
            !isempty(info.load_error) && (line *= "\tload_error=" * info.load_error)
            println(line)
        end
        return 0
    elseif sub == "register"
        return _cli_plugin_register(args)
    elseif sub == "remove"
        length(args) == 1 || throw(ArgumentError("`flopsy plugin remove` expects exactly one plugin name."))
        plugin_remove!(args[1])
        println("removed ", args[1])
        return 0
    end

    throw(ArgumentError("Unknown plugin subcommand `$sub`. Use `flopsy plugin --help`."))
end

function _cli_plugin_register(args::Vector{String})
    isempty(args) && throw(ArgumentError("`flopsy plugin register` expects a plugin name."))
    name = popfirst!(args)
    registry_url = nothing
    url = nothing
    path = nothing

    while !isempty(args)
        flag = popfirst!(args)
        if flag == "--registry"
            isempty(args) && throw(ArgumentError("`--registry` requires a URL value."))
            registry_url = popfirst!(args)
        elseif flag == "--url"
            isempty(args) && throw(ArgumentError("`--url` requires a package URL value."))
            url = popfirst!(args)
        elseif flag == "--path"
            isempty(args) && throw(ArgumentError("`--path` requires a local package path."))
            path = popfirst!(args)
        else
            throw(ArgumentError("Unknown plugin register flag `$flag`."))
        end
    end

    (url !== nothing && path !== nothing) &&
        throw(ArgumentError("Choose either `--url` or `--path`, not both."))

    plugin_register!(name; registry = registry_url, url = url, path = path)
    println("registered ", name)
    return 0
end

function _cli_xdmf(args::Vector{String})
    _wants_help(args) && return _cli_command_usage("xdmf")
    isempty(args) && throw(ArgumentError("`flopsy xdmf` expects an input HDF5 file."))

    input = popfirst!(args)
    output = nothing
    while !isempty(args)
        flag = popfirst!(args)
        flag == "--output" || throw(ArgumentError("Unknown xdmf flag `$flag`."))
        isempty(args) && throw(ArgumentError("`--output` requires a path value."))
        output = popfirst!(args)
    end
    out = write_xdmf_for_flopsy_h5(input, output)
    println(out)
    return 0
end

function _wants_help(args::Vector{String})
    return length(args) == 1 && args[1] in ("--help", "-h")
end

"""
    julia_main()::Cint

PackageCompiler-friendly entry point for the Flopsy CLI.
"""
function julia_main()::Cint
    return Cint(cli_main())
end
