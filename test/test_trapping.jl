@testset "Trapping models" begin

    @testset "FakeHotgatesModel rates" begin
        model = FakeHotgatesModel(2.0, 0.5)

        mobile  = [0.3]
        trapped = [0.2]
        defects = Float64[]
        dmobile  = zeros(1)
        dtrapped = zeros(1)

        hotgates_rates!(dmobile, dtrapped, model, mobile, defects, trapped, 300.0)

        trap_flux   = model.k_trap   * mobile[1] * (1 - trapped[1])
        detrap_flux = model.k_detrap * trapped[1]
        net = trap_flux - detrap_flux

        @test dmobile[1]  ≈ -net
        @test dtrapped[1] ≈  net

        # Conservation: total change should be zero
        @test dmobile[1] + dtrapped[1] ≈ 0.0 atol=1e-14
    end

    @testset "FakeHotgatesModel: length mismatch raises" begin
        model    = FakeHotgatesModel(1.0, 0.5)
        dmobile  = zeros(2)
        dtrapped = zeros(1)  # mismatch
        @test_throws ArgumentError hotgates_rates!(dmobile, dtrapped, model, [0.1, 0.2], Float64[], [0.1], 300.0)
    end

    @testset "SimpleTrappingReactionOperator: single node equilibrium" begin
        # At steady state: k_trap * c * (1-θ) = k_detrap * θ
        # → θ_eq = k_trap * c / (k_trap * c + k_detrap)
        k_trap   = 5.0
        k_detrap = 1.0
        c        = 0.4
        θ_eq     = k_trap * c / (k_trap * c + k_detrap)

        # Mesh1D requires nx ≥ 2; use two identical nodes
        nx     = 2
        mesh   = Mesh1D(1.0, nx)
        vars   = [VariableInfo(:c, :mobile, Set([:reaction])),
                  VariableInfo(:theta, :trap, Set([:reaction]))]
        layout = VariableLayout(vars)
        ctx    = SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())
        op     = SimpleTrappingReactionOperator(k_trap, k_detrap, 1, 2)

        # Both nodes at equilibrium
        u  = [c, θ_eq, c, θ_eq]   # layout: (var, node) column-major → [c1, θ1, c2, θ2]
        du = zeros(4)
        rhs!(du, op, u, ctx, 0.0)

        dU = state_view(du, layout, nx)
        @test dU[1, 1] ≈ 0.0 atol=1e-14   # dc/dt at node 1
        @test dU[2, 1] ≈ 0.0 atol=1e-14   # dθ/dt at node 1
    end

    @testset "SimpleTrappingReactionOperator: total conservation (no diffusion)" begin
        # d(c + θ)/dt = 0 when there is no diffusion
        k_trap   = 3.0
        k_detrap = 0.5
        nx       = 5
        mesh     = Mesh1D(1.0, nx)
        vars     = [VariableInfo(:c, :mobile, Set([:reaction])),
                    VariableInfo(:theta, :trap, Set([:reaction]))]
        layout   = VariableLayout(vars)
        ctx      = SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())
        op       = SimpleTrappingReactionOperator(k_trap, k_detrap, 1, 2)

        u  = rand(2 * nx)
        du = zeros(2 * nx)
        rhs!(du, op, u, ctx, 0.0)

        dU = state_view(du, layout, nx)
        # dc/dt + dθ/dt = 0 at every node
        for ix in 1:nx
            @test dU[1, ix] + dU[2, ix] ≈ 0.0 atol=1e-14
        end
    end

    @testset "build_trapping_model: end-to-end solve" begin
        mesh   = Mesh1D(1.0, 20)
        model  = build_trapping_model(
            mesh       = mesh,
            k_trap     = 2.0,
            k_detrap   = 0.5,
            diffusion_coefficient = 0.05,
        )

        nx     = length(mesh.x)
        nvars  = nvariables(model.layout)
        u0     = zeros(nvars * nx)
        U0     = state_view(u0, model.layout, nx)
        U0[1, nx ÷ 2] = 1.0   # pulse of mobile at centre

        cfg = SolverConfig(
            formulation = UnsplitFormulation(),
            algorithm   = Rodas5(autodiff=AutoFiniteDiff()),
            saveat      = range(0.0, 0.5; length=6),
        )

        sol    = solve_problem(model, u0, (0.0, 0.5), cfg)
        result = wrap_result(model, sol, Dict{String,Any}())

        @test sol.retcode == ReturnCode.Success
        @test length(result.solution.u) == 6

        # Total inventory (mobile + trap) with Neumann BCs should be conserved
        total_0 = sum(result.solution.u[1])
        total_f = sum(result.solution.u[end])
        @test total_0 ≈ total_f rtol=1e-4
    end

    @testset "build_trapping_model: analytic Jacobian solve" begin
        # Same model but the solver should use the analytic Jacobian automatically
        mesh  = Mesh1D(1.0, 15)
        model = build_trapping_model(
            mesh       = mesh,
            k_trap     = 1.0,
            k_detrap   = 0.2,
            diffusion_coefficient = 0.1,
        )
        @test supports_jacobian(OperatorSum(Tuple(active_operators(model))))

        nx    = length(mesh.x)
        nvars = nvariables(model.layout)
        u0    = zeros(nvars * nx)
        state_view(u0, model.layout, nx)[1, :] .= 0.5

        cfg = SolverConfig(
            formulation = UnsplitFormulation(),
            algorithm   = Rodas5(autodiff=AutoFiniteDiff()),
            saveat      = [0.0, 0.1],
        )
        sol = solve_problem(model, u0, (0.0, 0.1), cfg)
        @test sol.retcode == ReturnCode.Success
    end

    @testset "FakeHotgates via build_hotgates_trapping_model" begin
        nx     = 10
        mesh   = Mesh1D(1.0, nx)
        fake   = FakeHotgatesModel(3.0, 0.5)
        adaptor = HotgatesTrappingAdaptor(
            [1], [2],
            ["c"], ["theta"],
            String[],
            zeros(Float64, 0, nx),
        )

        model = build_hotgates_trapping_model(
            mesh                  = mesh,
            model                 = fake,
            adaptor               = adaptor,
            temperature           = ConstantTemperature(300.0),
            diffusion_coefficients = [0.05, 0.0],
        )

        nvars = nvariables(model.layout)
        u0    = zeros(nvars * nx)
        U0    = state_view(u0, model.layout, nx)
        U0[1, :] .= 0.3
        U0[2, :] .= 0.05

        cfg = SolverConfig(
            formulation = UnsplitFormulation(),
            algorithm   = Rodas5(autodiff=AutoFiniteDiff()),
            saveat      = [0.0, 0.2, 0.5],
        )
        sol = solve_problem(model, u0, (0.0, 0.5), cfg)
        @test sol.retcode == ReturnCode.Success

        # HotgatesReactionOperator doesn't support Jacobian — check this gracefully
        @test !supports_jacobian(OperatorSum(Tuple(active_operators(model))))
    end

    @testset "ToyReactionOperator: analytic decay" begin
        # u(t) = u0 * exp(-k*t), integrated over the domain
        k  = 0.5
        nx = 5
        mesh   = Mesh1D(1.0, nx)
        vars   = [VariableInfo(:u, :state, Set{Symbol}())]
        layout = VariableLayout(vars)
        react  = ToyReactionOperator(k)
        model  = build_rd_model(layout=layout, mesh=mesh, reaction=react)

        u0  = fill(1.0, nx)
        cfg = SolverConfig(
            formulation = UnsplitFormulation(),
            algorithm   = Rodas5(autodiff=AutoFiniteDiff()),
            saveat      = [0.0, 1.0, 2.0],
        )
        sol    = solve_problem(model, u0, (0.0, 2.0), cfg)
        result = wrap_result(model, sol, Dict{String,Any}())

        @test sol.retcode == ReturnCode.Success
        # Each node should follow exponential decay
        for ix in 1:nx
            @test sol.u[end][ix] ≈ exp(-k * 2.0) rtol=1e-4
        end
    end
end
