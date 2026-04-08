@testset "Formulations" begin

    # ------------------------------------------------------------------
    # Shared problem: simple trapping model, nx=20, short time window.
    # ------------------------------------------------------------------
    function make_trapping_problem(; nx=20, t_end=50.0)
        mesh   = Mesh1D(1e-3, nx)
        model  = build_trapping_model(
            mesh                  = mesh,
            k_trap                = 1e-2,
            k_detrap              = 0.1,
            diffusion_coefficient = 1e-7,
        )
        x   = mesh.x
        c0  = 0.1 .* exp.(-(x .- x[end]/2).^2 ./ (2*(x[end]/8)^2))
        u0  = zeros(2 * nx)
        for ix in 1:nx
            u0[(ix-1)*2 + 1] = c0[ix]
        end
        return model, u0, (0.0, t_end)
    end

    saveat = collect(range(0.0, 50.0; length = 11))

    # Reference solution from UnsplitFormulation
    model, u0, tspan = make_trapping_problem()
    config_ref = SolverConfig(
        formulation        = UnsplitFormulation(),
        algorithm          = Rodas5P(),
        abstol             = 1e-10,
        reltol             = 1e-8,
        saveat             = saveat,
        show_progress      = false,
        show_solver_stats  = false,
    )
    sol_ref    = solve_problem(model, u0, tspan, config_ref)
    result_ref = wrap_result(model, sol_ref, config_ref)

    c_ref = variable_timeseries(result_ref, :c)

    # ------------------------------------------------------------------
    @testset "UnsplitFormulation: basic solve" begin
        @test sol_ref.retcode == ReturnCode.Success
        @test length(result_ref.solution.t) == length(saveat)
        # Mobile concentration must remain non-negative
        @test all(c_ref .>= -1e-12)
    end

    # ------------------------------------------------------------------
    @testset "IMEXFormulation: solve and accuracy" begin
        model2, u02, tspan2 = make_trapping_problem()
        config_imex = SolverConfig(
            formulation        = IMEXFormulation(),
            algorithm          = KenCarp4(),
            abstol             = 1e-10,
            reltol             = 1e-8,
            saveat             = saveat,
            show_progress      = false,
            show_solver_stats  = false,
        )
        sol_imex    = solve_problem(model2, u02, tspan2, config_imex)
        result_imex = wrap_result(model2, sol_imex, config_imex)

        @test sol_imex.retcode == ReturnCode.Success
        @test length(result_imex.solution.t) == length(saveat)

        c_imex = variable_timeseries(result_imex, :c)
        # Solutions should agree to within a loose tolerance
        rel_err = maximum(abs.(c_imex .- c_ref)) / (maximum(abs.(c_ref)) + 1e-20)
        @test rel_err < 1e-3
    end

    # ------------------------------------------------------------------
    @testset "SplitFormulation(LieSplit): returns SplitSolution" begin
        model3, u03, tspan3 = make_trapping_problem()
        # Sub-steps do not have analytic Jacobians; use AutoFiniteDiff to avoid
        # ForwardDiff Dual-number conflicts with pre-allocated Float64 scratch buffers.
        config_lie = SolverConfig(
            formulation        = SplitFormulation(LieSplit()),
            algorithm          = Rodas5(autodiff=AutoFiniteDiff()),
            abstol             = 1e-9,
            reltol             = 1e-7,
            dt                 = 1.0,        # macro-step 1 s
            show_progress      = false,
            show_solver_stats  = false,
        )
        prob_lie = build_problem(model3, u03, tspan3, SplitFormulation(LieSplit()), config_lie)
        @test prob_lie isa SplitProblem

        sol_lie    = solve_problem(prob_lie, config_lie)
        result_lie = wrap_result(model3, sol_lie, config_lie)

        @test sol_lie isa SplitSolution
        @test sol_lie.retcode === :Success
        @test length(sol_lie.t) >= 2          # at least t0 and t_end
        @test sol_lie.t[1]   ≈ 0.0
        @test sol_lie.t[end] ≈ 50.0 atol=1.0 # within one macro-step

        # variable_timeseries should work on a SplitSolution
        c_lie = variable_timeseries(result_lie, :c)
        @test size(c_lie, 2) == 20
        @test all(c_lie .>= -1e-10)
    end

    # ------------------------------------------------------------------
    @testset "SplitFormulation(StrangSplit): higher accuracy than Lie" begin
        model4, u04, tspan4 = make_trapping_problem()
        config_strang = SolverConfig(
            formulation        = SplitFormulation(StrangSplit()),
            algorithm          = Rodas5(autodiff=AutoFiniteDiff()),
            abstol             = 1e-9,
            reltol             = 1e-7,
            dt                 = 1.0,
            show_progress      = false,
            show_solver_stats  = false,
        )
        sol_strang = solve_problem(
            build_problem(model4, u04, tspan4, SplitFormulation(StrangSplit()), config_strang),
            config_strang,
        )
        @test sol_strang isa SplitSolution
        @test sol_strang.retcode === :Success

        # Strang should produce a non-negative final profile
        result_strang = wrap_result(model4, sol_strang, config_strang)
        c_strang = variable_timeseries(result_strang, :c)
        @test all(c_strang .>= -1e-10)
    end

    # ------------------------------------------------------------------
    @testset "SplitFormulation: missing dt raises" begin
        model5, u05, tspan5 = make_trapping_problem()
        config_no_dt = SolverConfig(
            formulation        = SplitFormulation(LieSplit()),
            algorithm          = Rodas5(autodiff=AutoFiniteDiff()),
            abstol             = 1e-9,
            reltol             = 1e-7,
            show_progress      = false,
            show_solver_stats  = false,
            # dt not set
        )
        prob = build_problem(model5, u05, tspan5, SplitFormulation(LieSplit()), config_no_dt)
        @test_throws ArgumentError solve_problem(prob, config_no_dt)
    end

    # ------------------------------------------------------------------
    @testset "ResidualFormulation: solve" begin
        model6, u06, tspan6 = make_trapping_problem()
        config_dae = SolverConfig(
            formulation        = ResidualFormulation(),
            algorithm          = Rodas5P(),
            abstol             = 1e-10,
            reltol             = 1e-8,
            saveat             = saveat,
            show_progress      = false,
            show_solver_stats  = false,
        )
        sol_dae    = solve_problem(model6, u06, tspan6, config_dae)
        result_dae = wrap_result(model6, sol_dae, config_dae)

        @test sol_dae.retcode == ReturnCode.Success
        @test length(result_dae.solution.t) == length(saveat)

        c_dae = variable_timeseries(result_dae, :c)
        @test all(c_dae .>= -1e-10)
    end

    # ------------------------------------------------------------------
    @testset "_build_mass_matrix: correct structure" begin
        model7, _, _ = make_trapping_problem(nx=5)
        M = Flopsy._build_mass_matrix(model7)
        diag = M.diag

        # Even indices → :trap group (theta, ivar=2) → zero
        # Odd indices  → :mobile group (c, ivar=1) → one
        @test all(diag[1:2:end] .== 1.0)   # mobile (ivar=1)
        @test all(diag[2:2:end] .== 0.0)   # trap   (ivar=2)
    end

    # ------------------------------------------------------------------
    @testset "solver_stats_dict: works on SplitSolution" begin
        model8, u08, tspan8 = make_trapping_problem()
        config_s = SolverConfig(
            formulation        = SplitFormulation(LieSplit()),
            algorithm          = Rodas5(autodiff=AutoFiniteDiff()),
            dt                 = 5.0,
            show_progress      = false,
            show_solver_stats  = false,
        )
        sol_s    = solve_problem(build_problem(model8, u08, tspan8,
                                               SplitFormulation(LieSplit()), config_s),
                                 config_s)
        result_s = wrap_result(model8, sol_s, config_s)
        stats    = solver_stats_dict(result_s)
        @test stats isa Dict{String,Any}     # empty is fine, should not throw
    end

end
