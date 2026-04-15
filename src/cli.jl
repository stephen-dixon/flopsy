function cli_main(args::Vector{String} = copy(ARGS))
    isempty(args) && return _cli_usage(1)
    cmd = popfirst!(args)

    if cmd == "run"
        isempty(args) && return _cli_usage(1)
        run_input_deck(args[1])
        return 0
    elseif cmd == "validate"
        isempty(args) && return _cli_usage(1)
        validate_input_deck(args[1])
        println("valid")
        return 0
    elseif cmd == "syntax"
        return _cli_syntax(args)
    elseif cmd == "plugin"
        return _cli_plugin(args)
    elseif cmd == "xdmf"
        return _cli_xdmf(args)
    end

    return _cli_usage(1)
end

function _cli_usage(code::Int)
    println("usage: flopsy <run|validate|syntax|plugin|xdmf> ...")
    return code
end

function _cli_syntax(args::Vector{String})
    isempty(args) && return _cli_usage(1)
    sub = popfirst!(args)
    if sub == "list"
        for row in syntax_list()
            println("$(row[1])\t$(row[2])\tplugin=$(row[3])")
        end
        return 0
    elseif sub == "show"
        length(args) == 2 || return _cli_usage(1)
        info = syntax_show(Symbol(args[1]), Symbol(args[2]))
        println("domain: ", info.domain)
        println("type: ", info.type_name)
        println("plugin: ", info.plugin)
        println("help: ", info.help)
        println("parameters:")
        for param in info.parameters
            println("  - ", param.name, " required=", param.required, " default=", repr(param.default), " :: ", param.doc)
        end
        return 0
    end
    return _cli_usage(1)
end

function _cli_plugin(args::Vector{String})
    isempty(args) && return _cli_usage(1)
    sub = popfirst!(args)
    if sub == "list"
        for info in plugin_list()
            println(info.name, "\tsource=", info.source, isempty(info.load_error) ? "" : "\tload_error=$(info.load_error)")
        end
        return 0
    elseif sub == "register"
        isempty(args) && return _cli_usage(1)
        name = popfirst!(args)
        registry_url = nothing
        url = nothing
        path = nothing
        while !isempty(args)
            flag = popfirst!(args)
            if flag == "--registry"
                registry_url = popfirst!(args)
            elseif flag == "--url"
                url = popfirst!(args)
            elseif flag == "--path"
                path = popfirst!(args)
            else
                throw(ArgumentError("Unknown plugin register flag $flag"))
            end
        end
        plugin_register!(name; registry = registry_url, url = url, path = path)
        println("registered ", name)
        return 0
    elseif sub == "remove"
        isempty(args) && return _cli_usage(1)
        plugin_remove!(args[1])
        println("removed ", args[1])
        return 0
    end
    return _cli_usage(1)
end

function _cli_xdmf(args::Vector{String})
    isempty(args) && return _cli_usage(1)
    input = popfirst!(args)
    output = nothing
    while !isempty(args)
        flag = popfirst!(args)
        flag == "--output" || throw(ArgumentError("Unknown xdmf flag $flag"))
        output = popfirst!(args)
    end
    out = write_xdmf_for_flopsy_h5(input, output)
    println(out)
    return 0
end

function julia_main()::Cint
    try
        return Cint(cli_main())
    catch err
        println(stderr, sprint(showerror, err))
        return Cint(1)
    end
end
