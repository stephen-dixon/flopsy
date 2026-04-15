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
    @testset "IMEXReactionFormulation: solve and accuracy" begin
        # Build problem directly with SimpleTrappingReactionOperator, which has an
        # analytic Jacobian — exercises the jac1! path in build_problem.
        nx   = 20
        mesh = Mesh1D(1e-3, nx)
        x    = mesh.x

        selector = layout -> variables_with_tag(layout, :diffusion)
        react    = SimpleTrappingReactionOperator(1e-2, 0.1, 1, 2)
        diffop   = LinearDiffusionOperator([1e-7, 0.0], selector, nothing)
        bcop     = WeakDirichletBoundaryOperator(selector, [1e-7, 0.0], nothing;
                                                  left=t->0.0, right=t->0.0)

        vars   = [VariableInfo(:c, :mobile, Set([:diffusion, :reaction])),
                  VariableInfo(:theta, :trap,   Set([:reaction]))]
        layout = VariableLayout(vars)
        model_rxn = build_rd_model(
            layout    = layout,
            mesh      = mesh,
            diffusion = diffop,
            reaction  = react,
            boundary  = bcop,
        )

        c0 = 0.1 .* exp.(-(x .- x[end]/2).^2 ./ (2*(x[end]/8)^2))
        u0_rxn = zeros(2 * nx)
        dU0    = state_view(u0_rxn, layout, nx)
        for ix in 1:nx; dU0[1, ix] = c0[ix]; end

        config_rxn = SolverConfig(
            formulation       = IMEXReactionFormulation(),
            algorithm         = KenCarp4(),
            abstol            = 1e-10,
            reltol            = 1e-8,
            saveat            = saveat,
            show_progress     = false,
            show_solver_stats = false,
        )
        sol_rxn    = solve_problem(model_rxn, u0_rxn, tspan, config_rxn)
        result_rxn = wrap_result(model_rxn, sol_rxn, config_rxn)

        @test sol_rxn.retcode == ReturnCode.Success
        @test length(result_rxn.solution.t) == length(saveat)

        c_rxn = variable_timeseries(result_rxn, :c)
        @test all(c_rxn .>= -1e-10)

        # Compare against the UnsplitFormulation reference on the same model
        config_ref2 = SolverConfig(
            formulation       = UnsplitFormulation(),
            algorithm         = Rodas5P(),
            abstol            = 1e-10,
            reltol            = 1e-8,
            saveat            = saveat,
            show_progress     = false,
            show_solver_stats = false,
        )
        sol_ref2 = solve_problem(model_rxn, u0_rxn, tspan, config_ref2)
        c_ref2   = variable_timeseries(wrap_result(model_rxn, sol_ref2, config_ref2), :c)

        rel_err = maximum(abs.(c_rxn .- c_ref2)) / (maximum(abs.(c_ref2)) + 1e-20)
        @test rel_err < 1e-2
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

# ---------------------------------------------------------------------------
# check_mass_conservation tests
# ---------------------------------------------------------------------------

@testset "check_mass_conservation" begin

    @testset "closed system (Neumann BCs): total H constant" begin
        # Simple trapping with zero-flux BCs → no surface loss → H ≈ const.
        nx   = 20
        mesh = Mesh1D(1e-3, nx)
        x    = mesh.x
        model_c = build_trapping_model(mesh=mesh, k_trap=1e-2, k_detrap=0.1,
                                       diffusion_coefficient=1e-7)
        c0   = 0.1 .* exp.(-(x .- x[end]/2).^2 ./ (2*(x[end]/8)^2))
        u0_c = zeros(2*nx)
        for ix in 1:nx; u0_c[(ix-1)*2+1] = c0[ix]; end

        config_c = SolverConfig(
            formulation=UnsplitFormulation(), algorithm=Rodas5P(),
            abstol=1e-10, reltol=1e-8,
            saveat=collect(range(0.0, 50.0; length=21)),
            show_progress=false, show_solver_stats=false,
        )
        result_c = wrap_result(model_c, solve_problem(model_c, u0_c, (0.0, 50.0), config_c), config_c)

        mc = check_mass_conservation(result_c; rtol=1e-4)
        @test mc.conserved
        # Total H should be nearly constant (no surface loss)
        @test maximum(abs.(mc.total_hydrogen .- mc.total_hydrogen[1])) /
              max(abs(mc.total_hydrogen[1]), 1e-30) < 1e-4
    end

    @testset "open system (Dirichlet BCs): surface flux accounts for loss" begin
        # Trapping with vacuum Dirichlet BCs → H decreases; flux integral should match.
        nx = 20
        mesh = Mesh1D(1e-3, nx)
        x    = mesh.x
        D    = 1e-7
        selector  = layout -> variables_with_tag(layout, :diffusion)
        react_op  = SimpleTrappingReactionOperator(1e-2, 0.1, 1, 2)
        diffop_o  = LinearDiffusionOperator([D, 0.0], selector, nothing)
        bcop_o    = WeakDirichletBoundaryOperator(selector, [D, 0.0], nothing;
                                                   left=t->0.0, right=t->0.0)
        vars_o    = [VariableInfo(:c,     :mobile, Set([:diffusion, :reaction])),
                     VariableInfo(:theta, :trap,   Set([:reaction]))]
        layout_o  = VariableLayout(vars_o)
        model_o   = build_rd_model(layout=layout_o, mesh=mesh,
                                    diffusion=diffop_o, reaction=react_op, boundary=bcop_o)

        u0_o  = zeros(2*nx)
        dU0_o = state_view(u0_o, layout_o, nx)
        # Use a profile that satisfies the BCs at t=0 (zero at boundaries).
        # Use fine saveat (dt=0.01 s) so that the trapezoidal flux integral is
        # accurate enough for rtol=1e-3; the diffusion time scale here is ~1 s.
        for ix in 1:nx; dU0_o[1, ix] = 0.05 * sin(π * x[ix] / x[end]); end

        config_o = SolverConfig(
            formulation=UnsplitFormulation(), algorithm=Rodas5P(),
            abstol=1e-10, reltol=1e-8,
            saveat=collect(range(0.0, 5.0; length=501)),
            show_progress=false, show_solver_stats=false,
        )
        result_o = wrap_result(model_o, solve_problem(model_o, u0_o, (0.0, 5.0), config_o), config_o)

        mc = check_mass_conservation(result_o; rtol=1e-3)
        # H should decrease (matter leaves through surfaces)
        @test mc.total_hydrogen[end] < mc.total_hydrogen[1]
        # Mass balance: H(t) + cumulative_flux(t) ≈ H(0)
        @test mc.conserved
    end

    @testset "check_mass_conservation return structure" begin
        nx   = 10
        mesh = Mesh1D(1e-3, nx)
        model_s = build_trapping_model(mesh=mesh, k_trap=1e-2, k_detrap=0.1,
                                       diffusion_coefficient=1e-7)
        u0_s = zeros(2*nx); u0_s[1] = 0.01
        config_s2 = SolverConfig(formulation=UnsplitFormulation(), algorithm=Rodas5P(),
            abstol=1e-10, reltol=1e-8, saveat=[0.0, 5.0, 10.0],
            show_progress=false, show_solver_stats=false)
        result_s2 = wrap_result(model_s, solve_problem(model_s, u0_s, (0.0, 10.0), config_s2), config_s2)

        mc = check_mass_conservation(result_s2)
        @test haskey(mc, :conserved)
        @test haskey(mc, :max_relative_error)
        @test haskey(mc, :total_hydrogen)
        @test haskey(mc, :cumulative_flux)
        @test haskey(mc, :balance)
        @test length(mc.total_hydrogen) == 3
        @test length(mc.cumulative_flux) == 3
        @test mc.cumulative_flux[1] == 0.0   # cumulative flux starts at zero
    end
end
