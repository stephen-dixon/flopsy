const _CLI_MAIN_HELP = """
flopsy: registry-driven reaction-diffusion input-deck runner

Usage:
  flopsy run <input.toml> [--problem <name>] [--output-dir <dir>]
  flopsy validate <input.toml> [--problem <name>]
  flopsy report <fields.h5> [--csv <path>]
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
  report    Print key metrics from an HDF5 output file.
  syntax    Inspect built-in and plugin-provided registered syntax.
  plugin    Manage runtime plugins loaded into the syntax registry.
  xdmf      Generate an XDMF companion for an HDF5 field output file.

Use `flopsy <command> --help` for command-specific help.
"""

const _CLI_COMMAND_HELP = Dict(
    "run" => """
Usage:
  flopsy run <input.toml> [--problem <name>] [--output-dir <dir>]

Run a registry-driven input deck from a TOML file.

Options:
  --problem <name>     Select a specific [problem.<name>] block (default: alphabetically first).
  --output-dir <dir>   Write all outputs under this directory instead of their declared paths.
""",
    "validate" => """
Usage:
  flopsy validate <input.toml> [--problem <name>]

Validate a registry-driven input deck without running the solver.

Options:
  --problem <name>  Select a specific [problem.<name>] block to validate (default: alphabetically first).
""",
    "report" => """
Usage:
  flopsy report <fields.h5> [--csv <path>]

Print key metrics from a Flopsy HDF5 field output file.
If a companion summary CSV with the same base name exists, it is loaded for flux data.

Options:
  --csv <path>  Path to write a summary CSV if one is not already present.
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
"""
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
    elseif cmd == "report"
        return _cli_report(args)
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

    path = nothing
    problem_name = nothing
    output_dir = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--problem"
            i += 1
            i > length(args) && throw(ArgumentError("`--problem` requires a block name."))
            problem_name = args[i]
        elseif arg == "--output-dir"
            i += 1
            i > length(args) && throw(ArgumentError("`--output-dir` requires a directory path."))
            output_dir = args[i]
        elseif startswith(arg, "-")
            throw(ArgumentError("Unknown flag `$arg`. Use `flopsy run --help` for usage."))
        elseif path === nothing
            path = arg
        else
            throw(ArgumentError("`flopsy run` accepts only one input deck path."))
        end
        i += 1
    end

    path === nothing && throw(ArgumentError("`flopsy run` requires an input deck path."))
    run_input_deck(path; problem_name = problem_name, output_dir = output_dir)
    return 0
end

function _cli_validate(args::Vector{String})
    _wants_help(args) && return _cli_command_usage("validate")

    path = nothing
    problem_name = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--problem"
            i += 1
            i > length(args) && throw(ArgumentError("`--problem` requires a block name."))
            problem_name = args[i]
        elseif startswith(arg, "-")
            throw(ArgumentError("Unknown flag `$arg`. Use `flopsy validate --help` for usage."))
        elseif path === nothing
            path = arg
        else
            throw(ArgumentError("`flopsy validate` accepts only one input deck path."))
        end
        i += 1
    end

    path === nothing && throw(ArgumentError("`flopsy validate` requires an input deck path."))
    validate_input_deck(path; problem_name = problem_name)
    println("valid")
    return 0
end

function _cli_report(args::Vector{String})
    _wants_help(args) && return _cli_command_usage("report")
    isempty(args) && throw(ArgumentError("`flopsy report` expects an HDF5 file path."))

    h5_path = popfirst!(args)
    csv_path = nothing

    while !isempty(args)
        flag = popfirst!(args)
        if flag == "--csv"
            isempty(args) && throw(ArgumentError("`--csv` requires a path value."))
            csv_path = popfirst!(args)
        else
            throw(ArgumentError("Unknown report flag `$flag`. Use `flopsy report --help`."))
        end
    end

    _report_h5(h5_path; csv_path = csv_path)
    return 0
end

function _report_h5(h5_path::AbstractString; csv_path = nothing)
    isfile(h5_path) || throw(ArgumentError("File not found: $h5_path"))

    h5open(h5_path, "r") do h5
        t = read(h5["time"])
        nt = length(t)
        fields = sort!(collect(String.(keys(h5["fields"]))))
        x = read(h5["mesh/x"])
        dx = read(h5["mesh/dx"])
        nx = length(x)

        println("HDF5 file   : ", h5_path)
        println("Time steps  : ", nt)
        println("Time range  : [", t[1], ", ", t[end], "]")
        println("Nodes       : ", nx)
        println("Fields      : ", join(fields, ", "))
        println()
        println("Summary at final time (t = ", t[end], "):")

        for fname in fields
            data = read(h5["fields"][fname])
            final = @view data[end, :]
            total = sum(final) * dx
            peak = maximum(abs, final)
            @printf("  %-24s  total = %10.4e  peak = %10.4e\n", fname, total, peak)
        end
    end

    # Check for a companion summary CSV
    auto_csv = replace(h5_path, r"\.h5$" => "_summary.csv")
    if isfile(auto_csv) && csv_path === nothing
        println()
        println("Companion CSV: ", auto_csv, " (use `flopsy report --csv <path>` to specify)")
    end

    if csv_path !== nothing && !isfile(csv_path)
        println()
        println("Note: CSV generation requires a SimulationResult. Re-run with a [output.*.type = \"summary_csv\"] block.")
    end

    return nothing
end

function _cli_syntax(args::Vector{String})
    isempty(args) && return _cli_command_usage("syntax")
    _wants_help(args) && return _cli_command_usage("syntax")

    sub = popfirst!(args)
    if sub == "list"
        isempty(args) ||
            throw(ArgumentError("`flopsy syntax list` does not accept extra arguments."))
        for row in syntax_list()
            println("$(row[1])\t$(row[2])\tplugin=$(row[3])")
        end
        return 0
    elseif sub == "show"
        length(args) == 2 ||
            throw(ArgumentError("`flopsy syntax show` expects <domain> <type>."))
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
            param.allowed_values === nothing ||
                print(" allowed=", repr(param.allowed_values))
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
        isempty(args) ||
            throw(ArgumentError("`flopsy plugin list` does not accept extra arguments."))
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
        length(args) == 1 ||
            throw(ArgumentError("`flopsy plugin remove` expects exactly one plugin name."))
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
