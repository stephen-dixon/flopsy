@testset "Config and Factory" begin
    @testset "parse legacy config into typed config" begin
        raw = Dict(
            "model_type" => "toy_trapping",
            "length" => 1.0,
            "nx" => 21,
            "tspan" => [0.0, 2.0],
            "k_trap" => 2.0,
            "k_detrap" => 0.25,
            "diffusion_coefficient" => 0.05,
            "saveat" => [0.0, 0.5, 2.0],
            "algorithm" => "Rodas5"
        )

        cfg = parse_config(raw)

        @test cfg isa ProblemConfig
        @test cfg.problem_type === :trapping_1d
        @test cfg.mesh.nx == 21
        @test cfg.solver.algorithm === :Rodas5
        @test cfg.solver.saveat == [0.0, 0.5, 2.0]
        @test getproperty(cfg.parameters, :k_trap) == 2.0
    end

    @testset "validation fails early on bad split config" begin
        cfg = ProblemConfig(
            :diffusion_1d,
            MeshConfig(:uniform_1d, 11, 0.0, 1.0),
            InputSolverConfig(:split, :Rodas5, nothing, 1e-8, 1e-6, [0.0, 1.0]),
            BoundaryConditionConfig{Float64}[],
            (; t0 = 0.0, tend = 1.0, diffusion_coefficient = 0.1)
        )

        @test_throws ArgumentError validate(cfg)
    end

    @testset "problem factory builds executable trapping problem" begin
        cfg = ProblemConfig(
            :trapping_1d,
            MeshConfig(:uniform_1d, 17, 0.0, 1.0),
            InputSolverConfig(:unsplit, :Rodas5, nothing, 1e-8, 1e-6, [0.0, 0.1]),
            BoundaryConditionConfig{Float64}[],
            (; t0 = 0.0, tend = 0.1, k_trap = 2.0, k_detrap = 0.25,
                diffusion_coefficient = 0.05, initial_mobile_pulse_amplitude = 1.0,
                initial_trap_occupancy = 0.0)
        )

        problem = build_problem(cfg)
        result = solve(problem)

        @test problem isa SimulationProblem
        @test result.solution.t == [0.0, 0.1]
        @test result.metadata["problem_type"] == "trapping_1d"
    end

    @testset "load_config reads typed TOML deck" begin
        cfg = load_config(joinpath(@__DIR__, "..", "examples", "trapping_1d.toml"))
        @test cfg.problem_type === :trapping_1d
        @test cfg.mesh.kind === :uniform_1d
        @test cfg.solver.formulation === :unsplit
    end
end
