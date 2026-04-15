@testset "Diffusion operators" begin

    # Helpers
    function make_ctx(nx; nvars = 1)
        mesh = Mesh1D(1.0, nx)
        vars = [VariableInfo(Symbol("u$i"), :state, Set([:diffusion])) for i in 1:nvars]
        layout = VariableLayout(vars)
        SystemContext(layout, nx, mesh, Dict{Symbol, Any}(), Dict{Symbol, Any}())
    end

    selector_all = layout -> collect(1:nvariables(layout))

    @testset "LinearDiffusionOperator: flat profile → zero flux" begin
        nx = 10
        ctx = make_ctx(nx)
        D = 1.0
        op = LinearDiffusionOperator([D], selector_all, nothing)

        u = ones(nx)
        du = zeros(nx)
        rhs!(du, op, u, ctx, 0.0)
        @test all(≈(0.0; atol = 1e-14), du)
    end

    @testset "LinearDiffusionOperator: interior stencil value" begin
        # Quadratic profile u = x^2 → d²u/dx² = 2 everywhere
        nx = 11
        ctx = make_ctx(nx)
        dx = ctx.mesh.dx
        x = ctx.mesh.x
        D = 1.0
        op = LinearDiffusionOperator([D], selector_all, nothing)

        u = x .^ 2
        du = zeros(nx)
        rhs!(du, op, u, ctx, 0.0)

        # Interior nodes should each give ≈ D * 2 = 2.0
        @test all(≈(2.0; atol = 1e-10), du[2:(nx - 1)])
    end

    @testset "DirichletBoundaryOperator: boundary correction" begin
        nx = 5
        ctx = make_ctx(nx)
        dx = ctx.mesh.dx
        D = 2.0
        g_left = 1.0
        g_right = 3.0

        op = WeakDirichletBoundaryOperator(selector_all, [D], nothing;
            left = t -> g_left, right = t -> g_right)

        u = zeros(nx)
        du = zeros(nx)
        rhs!(du, op, u, ctx, 0.0)

        # Only boundary nodes should be non-zero
        @test du[1] ≈ D * (g_left - u[1]) / dx^2
        @test du[nx] ≈ D * (g_right - u[nx]) / dx^2
        @test all(≈(0.0; atol = 1e-14), du[2:(nx - 1)])
    end

    @testset "DirichletBoundaryOperator: only left BC" begin
        nx = 5
        ctx = make_ctx(nx)
        dx = ctx.mesh.dx
        D = 1.0

        op = WeakDirichletBoundaryOperator(selector_all, [D], nothing; left = t -> 1.0)

        u = zeros(nx)
        du = zeros(nx)
        rhs!(du, op, u, ctx, 0.0)

        @test du[1] ≈ D / dx^2
        @test du[nx] ≈ 0.0
    end

    @testset "Combined Neumann + Dirichlet gives centred stencil at boundary" begin
        nx = 5
        ctx = make_ctx(nx)
        dx = ctx.mesh.dx
        D = 1.0
        g = 0.5

        diffop = LinearDiffusionOperator([D], selector_all, nothing)
        boundop = WeakDirichletBoundaryOperator(selector_all, [D], nothing; left = t -> g)

        u = collect(Float64, 1:nx)
        du = zeros(nx)
        # Apply both operators (accumulating)
        rhs!(du, diffop, u, ctx, 0.0)
        rhs!(du, boundop, u, ctx, 0.0)

        # Node 1: Neumann gives D*(u[2]-u[1])/dx²; Dirichlet adds D*(g-u[1])/dx²
        # Combined: D*(u[2] - 2*u[1] + g)/dx² — standard centred Dirichlet stencil
        expected = D * (u[2] - 2*u[1] + g) / dx^2
        @test du[1] ≈ expected atol=1e-12
    end

    @testset "NullOperator Jacobian" begin
        nx = 3
        ctx = make_ctx(nx)
        n = nx
        J = zeros(n, n)
        op = NullOperator()
        u = ones(n)
        jacobian!(J, op, u, ctx, 0.0)
        @test all(J .== 0.0)  # NullOperator leaves J unchanged (zero)
    end

    @testset "supports_jacobian flags" begin
        selector = layout -> [1]
        @test supports_jacobian(LinearDiffusionOperator([1.0], selector, nothing))
        @test supports_jacobian(WeakDirichletBoundaryOperator(selector, [1.0], nothing; left = t->0.0))
        @test supports_jacobian(NullOperator())
        @test supports_jacobian(ToyReactionOperator(0.1))
        @test supports_jacobian(SimpleTrappingReactionOperator(1.0, 0.5, 1, 2))
    end
end

@testset "HotgatesReactionOperator does not support Jacobian" begin
    fake = FakeHotgatesModel(1.0, 0.5)
    nx = 3
    adaptor = HotgatesTrappingAdaptor([1], [2], ["c"], ["theta"], String[], zeros(0, nx))
    op = HotgatesReactionOperator(fake, adaptor, ConstantTemperature(300.0))
    @test !supports_jacobian(op)
end

@testset "OperatorSum Jacobian flag" begin
    selector = layout -> variables_with_tag(layout, :diffusion)

    @testset "all support → OperatorSum supports" begin
        diffop = LinearDiffusionOperator([1.0], selector, nothing)
        react = ToyReactionOperator(0.1)
        total = OperatorSum((diffop, react, NullOperator()))
        @test supports_jacobian(total)
    end

    @testset "one missing → OperatorSum does not support" begin
        diffop = LinearDiffusionOperator([1.0], selector, nothing)
        fake = FakeHotgatesModel(1.0, 0.5)
        nx = 3
        adaptor = HotgatesTrappingAdaptor(
            [1], [2], ["c"], ["theta"], String[], zeros(0, nx))
        hotop = HotgatesReactionOperator(fake, adaptor, ConstantTemperature(300.0))
        total = OperatorSum((diffop, hotop))
        @test !supports_jacobian(total)
    end
end
