@testset "Registry and CLI" begin
    @testset "built-in syntax is registered through the registry" begin
        registry = build_registry()
        mesh_spec = lookup_syntax(registry, :mesh, :uniform_1d)
        backend_spec = lookup_syntax(registry, :backend, :trapping_1d)
        bc_spec = lookup_syntax(registry, :bc, :dirichlet)

        @test mesh_spec.domain === :mesh
        @test backend_spec.domain === :backend
        @test bc_spec.type_name === :dirichlet
        @test_throws ArgumentError lookup_syntax(registry, :backend, :does_not_exist)
    end

    @testset "named-block input deck builds and solves" begin
        result = run_input_deck(joinpath(@__DIR__, "..", "examples", "trapping_1d_with_bc_ic.toml"))
        @test result.solution.t == [0.0, 0.25, 0.5, 0.75, 1.0]
    end

    @testset "omitted BCs give closed zero-flux behaviour" begin
        result = run_input_deck(joinpath(@__DIR__, "..", "examples", "closed_system_no_bc.toml"))
        inventory = integrated_variable(result, :u)
        @test all(isapprox.(inventory, fill(inventory[1], length(inventory)); rtol = 1e-6, atol = 1e-8))
    end

    @testset "species-targeted IC overlap is rejected" begin
        raw = Dict(
            "mesh" => Dict("main" => Dict("type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 11)),
            "backend" => Dict("main" => Dict("type" => "diffusion_1d", "diffusion_coefficient" => 0.1)),
            "ic" => Dict(
                "a" => Dict("type" => "uniform_species", "species" => "u", "value" => 1.0),
                "b" => Dict("type" => "uniform_species", "species" => "u", "value" => 2.0),
            ),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "main",
                "backend" => "main",
                "ics" => ["a", "b"],
                "tspan" => [0.0, 1.0],
            )),
        )
        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        @test_throws ArgumentError Flopsy.build_simulation(ctx)
    end

    @testset "species-targeted BC validation rejects stationary species" begin
        raw = Dict(
            "mesh" => Dict("main" => Dict("type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 11)),
            "backend" => Dict("main" => Dict(
                "type" => "trapping_1d",
                "k_trap" => 1.0,
                "k_detrap" => 0.1,
                "diffusion_coefficient" => 0.01,
            )),
            "bc" => Dict("trap" => Dict(
                "type" => "dirichlet",
                "species" => "H_trapped",
                "boundary" => "left",
                "value" => 0.0,
            )),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "main",
                "backend" => "main",
                "bcs" => ["trap"],
                "tspan" => [0.0, 1.0],
            )),
        )
        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        @test_throws ArgumentError Flopsy.build_simulation(ctx)
    end

    @testset "CLI validate syntax xdmf and run" begin
        deck = joinpath(@__DIR__, "..", "examples", "trapping_1d_with_bc_ic.toml")
        tmp = mktempdir()
        cd(tmp) do
            @test cli_main(["validate", deck]) == 0
            @test cli_main(["syntax", "list"]) == 0
            @test cli_main(["syntax", "show", "bc", "dirichlet"]) == 0
            @test cli_main(["run", deck]) == 0
            h5 = joinpath(tmp, "trapping_fields.h5")
            @test isfile(h5)
            @test cli_main(["xdmf", h5, "--output", joinpath(tmp, "manual.xdmf")]) == 0
            @test isfile(joinpath(tmp, "manual.xdmf"))
            @test isfile(joinpath(tmp, "trapping_fields.xdmf"))
        end
    end

    @testset "plugin install list and remove use managed environment" begin
        tmp_home = mktempdir()
        plugin_src = mktempdir()
        name = "DummyFlopsyPlugin"
        mkpath(joinpath(plugin_src, "src"))
        write(joinpath(plugin_src, "Project.toml"), """
name = \"$name\"
uuid = \"1f74d446-c9be-4a9b-a650-79a6fa7a8f11\"
version = \"0.1.0\"
""")
        write(joinpath(plugin_src, "src", "$name.jl"), """
module $name
function register_flopsy_plugin!(registry)
    return registry
end
end
""")

        withenv("HOME" => tmp_home) do
            @test isempty(plugin_list())
            plugin_register!(name; path = plugin_src)
            infos = plugin_list()
            @test length(infos) == 1
            @test infos[1].name == name
            registry = build_registry()
            @test registry isa SyntaxRegistry
            plugin_remove!(name)
            @test isempty(plugin_list())
        end
    end
end
