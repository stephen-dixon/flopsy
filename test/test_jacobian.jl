@testset "Analytic Jacobians" begin

    # Finite-difference Jacobian for a standalone operator.
    # f!(du, u, ctx) must zero du then accumulate, or be tested via OperatorSum.
    function fd_jacobian(f!, u0, ctx)
        n  = length(u0)
        J  = zeros(n, n)
        du0 = zeros(n)
        f!(du0, copy(u0), ctx, 0.0)
        ε = sqrt(eps(Float64))
        for j in 1:n
            u_p = copy(u0)
            u_p[j] += ε
            du_p = zeros(n)
            f!(du_p, u_p, ctx, 0.0)
            J[:, j] .= (du_p .- du0) ./ ε
        end
        return J
    end

    # Wrap single operator into a zero-initialising f! (like OperatorSum does)
    function op_f!(op)
        return (du, u, ctx, t) -> begin
            fill!(du, 0.0)
            rhs!(du, op, u, ctx, t)
        end
    end

    function make_ctx(nx; nvars=1)
        mesh = Mesh1D(1.0, nx)
        vars = [VariableInfo(Symbol("u$i"), :state, Set([:diffusion])) for i in 1:nvars]
        layout = VariableLayout(vars)
        SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())
    end

    selector_all = layout -> collect(1:nvariables(layout))

    @testset "LinearDiffusionOperator — single variable" begin
        nx  = 6
        ctx = make_ctx(nx)
        D   = 2.5
        op  = LinearDiffusionOperator([D], selector_all, nothing)

        u0 = randn(nx)
        J_fd = fd_jacobian(op_f!(op), u0, ctx)

        J_an = zeros(nx, nx)
        jacobian!(J_an, op, u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "LinearDiffusionOperator — two variables, one diffuses" begin
        nx    = 5
        nvars = 2
        ctx   = make_ctx(nx; nvars=nvars)
        n     = nvars * nx
        # only variable 1 diffuses
        sel   = layout -> [1]
        op    = LinearDiffusionOperator([1.5, 0.0], sel, nothing)

        u0 = randn(n)
        J_fd = fd_jacobian(op_f!(op), u0, ctx)
        J_an = zeros(n, n)
        jacobian!(J_an, op, u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "WeakDirichletBoundaryOperator — both BCs" begin
        nx  = 6
        ctx = make_ctx(nx)
        D   = 1.0
        op  = WeakDirichletBoundaryOperator(selector_all, [D], nothing;
                                             left=t->1.0, right=t->2.0)

        u0 = randn(nx)
        J_fd = fd_jacobian(op_f!(op), u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, op, u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "ToyReactionOperator" begin
        nx  = 5
        ctx = make_ctx(nx)
        op  = ToyReactionOperator(0.3)

        u0 = randn(nx)
        J_fd = fd_jacobian(op_f!(op), u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, op, u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "SimpleTrappingReactionOperator" begin
        nx    = 4
        nvars = 2
        ctx   = make_ctx(nx; nvars=nvars)
        n     = nvars * nx
        op    = SimpleTrappingReactionOperator(2.0, 0.5, 1, 2)

        u0 = abs.(randn(n)) .* 0.4  # keep mobile/trap in (0,1) range
        u0 = clamp.(u0, 0.01, 0.99)

        J_fd = fd_jacobian(op_f!(op), u0, ctx)
        J_an = zeros(n, n)
        jacobian!(J_an, op, u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "OperatorSum Jacobian (diffusion + reaction)" begin
        nx    = 5
        nvars = 2
        ctx   = make_ctx(nx; nvars=nvars)
        n     = nvars * nx

        sel  = layout -> [1]  # only variable 1 diffuses
        diff = LinearDiffusionOperator([1.2, 0.0], sel, nothing)
        react = SimpleTrappingReactionOperator(3.0, 0.8, 1, 2)
        total = OperatorSum((diff, react))

        @test supports_jacobian(total)

        u0 = clamp.(abs.(randn(n)) .* 0.5, 0.01, 0.99)
        u0[1:2:end] .= clamp.(u0[1:2:end], 0.01, 5.0)   # mobile can exceed 1

        f_sum! = (du, u, ctx, t) -> (fill!(du, 0.0); rhs!(du, total, u, ctx, t))
        J_fd = fd_jacobian(f_sum!, u0, ctx)

        J_an = zeros(n, n)
        jacobian!(J_an, total, u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-4
    end

    @testset "Sparsity prototype has correct structure" begin
        nx    = 4
        nvars = 2
        ctx   = make_ctx(nx; nvars=nvars)
        n     = nvars * nx

        sel    = layout -> [1]
        diff   = LinearDiffusionOperator([1.0, 0.0], sel, nothing)
        react  = SimpleTrappingReactionOperator(1.0, 0.5, 1, 2)
        ops    = (diff, react)
        layout = ctx.layout

        proto = Flopsy._build_jac_prototype(
            SystemModel(layout, (;), ctx),
            ops,
        )

        @test size(proto) == (n, n)

        # Analytic Jacobian entries must be a subset of the nonzero pattern
        u0 = clamp.(abs.(randn(n)) .* 0.5, 0.01, 0.99)
        J_an = zeros(n, n)
        jacobian!(J_an, OperatorSum(ops), u0, ctx, 0.0)

        # Every entry of the analytic Jacobian must correspond to a nonzero in the prototype
        I_proto, J_proto, _ = findnz(proto)
        proto_set = Set(zip(I_proto, J_proto))
        for i in 1:n, j in 1:n
            if J_an[i, j] != 0.0
                @test (i, j) ∈ proto_set
            end
        end
    end
end

@testset "Strong DirichletBoundaryOperator Jacobians" begin

    function fd_jacobian(f!, u0, ctx)
        n   = length(u0)
        J   = zeros(n, n)
        du0 = zeros(n)
        f!(du0, copy(u0), ctx, 0.0)
        ε = sqrt(eps(Float64))
        for j in 1:n
            u_p  = copy(u0)
            u_p[j] += ε
            du_p = zeros(n)
            f!(du_p, u_p, ctx, 0.0)
            J[:, j] .= (du_p .- du0) ./ ε
        end
        return J
    end

    function make_ctx(nx; nvars=1)
        mesh = Mesh1D(1.0, nx)
        vars = [VariableInfo(Symbol("u$i"), :state, Set([:diffusion])) for i in 1:nvars]
        layout = VariableLayout(vars)
        SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())
    end

    # Combined diffusion + strong Dirichlet f! helper
    function combined_f!(diffop, bcop)
        return (du, u, ctx, t) -> begin
            fill!(du, 0.0)
            rhs!(du, diffop, u, ctx, t)
            rhs!(du, bcop,   u, ctx, t)
        end
    end

    nx       = 6
    ctx      = make_ctx(nx)
    D        = 1.5
    selector = layout -> collect(1:nvariables(layout))
    diffop   = LinearDiffusionOperator([D], selector, nothing)
    coeffs   = [D]

    @testset "PenaltyMethod — left BC" begin
        op = DirichletBoundaryOperator(:left, t -> 1.0, selector; method = PenaltyMethod(1e6))
        f! = combined_f!(diffop, op)

        u0   = randn(nx)
        J_fd = fd_jacobian(f!, u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, diffop, u0, ctx, 0.0)
        jacobian!(J_an, op,     u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "PenaltyMethod — right BC" begin
        op = DirichletBoundaryOperator(:right, t -> 2.0, selector; method = PenaltyMethod(1e6))
        f! = combined_f!(diffop, op)

        u0   = randn(nx)
        J_fd = fd_jacobian(f!, u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, diffop, u0, ctx, 0.0)
        jacobian!(J_an, op,     u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "MassMatrixMethod — left BC" begin
        op = DirichletBoundaryOperator(:left, t -> 0.5, selector; method = MassMatrixMethod())
        f! = combined_f!(diffop, op)

        u0   = randn(nx)
        J_fd = fd_jacobian(f!, u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, diffop, u0, ctx, 0.0)
        jacobian!(J_an, op,     u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "MassMatrixMethod — right BC" begin
        op = DirichletBoundaryOperator(:right, t -> 0.0, selector; method = MassMatrixMethod())
        f! = combined_f!(diffop, op)

        u0   = randn(nx)
        J_fd = fd_jacobian(f!, u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, diffop, u0, ctx, 0.0)
        jacobian!(J_an, op,     u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "CallbackMethod — no Jacobian contribution" begin
        # CallbackMethod rhs! is a no-op; its Jacobian should be zero.
        op = DirichletBoundaryOperator(:left, t -> 1.0, selector; method = CallbackMethod())

        u0   = randn(nx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, op, u0, ctx, 0.0)

        @test all(J_an .== 0.0)
    end

    @testset "EliminatedMethod — left BC (combined with diffusion)" begin
        op = DirichletBoundaryOperator(:left, t -> 0.0, selector, coeffs; method = EliminatedMethod())
        f! = combined_f!(diffop, op)

        u0   = randn(nx)
        J_fd = fd_jacobian(f!, u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, diffop, u0, ctx, 0.0)
        jacobian!(J_an, op,     u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "EliminatedMethod — right BC (combined with diffusion)" begin
        op = DirichletBoundaryOperator(:right, t -> 0.0, selector, coeffs; method = EliminatedMethod())
        f! = combined_f!(diffop, op)

        u0   = randn(nx)
        J_fd = fd_jacobian(f!, u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, diffop, u0, ctx, 0.0)
        jacobian!(J_an, op,     u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end

    @testset "OperatorSum — diffusion + PenaltyMethod both sides" begin
        left_bc  = DirichletBoundaryOperator(:left,  t -> 1.0, selector; method = PenaltyMethod(1e5))
        right_bc = DirichletBoundaryOperator(:right, t -> 0.0, selector; method = PenaltyMethod(1e5))
        total    = OperatorSum((diffop, left_bc, right_bc))

        @test supports_jacobian(total)

        u0   = randn(nx)
        f!   = (du, u, ctx, t) -> (fill!(du, 0.0); rhs!(du, total, u, ctx, t))
        J_fd = fd_jacobian(f!, u0, ctx)
        J_an = zeros(nx, nx)
        jacobian!(J_an, total, u0, ctx, 0.0)

        @test norm(J_an - J_fd) / (norm(J_fd) + 1e-12) < 1e-5
    end
end

@testset "FakeHotgatesModel finite-difference Jacobian fallback" begin

    function fd_jacobian(f!, u0, ctx)
        n   = length(u0)
        J   = zeros(n, n)
        du0 = zeros(n)
        f!(du0, copy(u0), ctx, 0.0)
        ε = sqrt(eps(Float64))
        for j in 1:n
            u_p  = copy(u0)
            u_p[j] += ε
            du_p = zeros(n)
            f!(du_p, u_p, ctx, 0.0)
            J[:, j] .= (du_p .- du0) ./ ε
        end
        return J
    end

    function make_ctx(nx; nvars=1)
        mesh = Mesh1D(1.0, nx)
        vars = [VariableInfo(Symbol("u$i"), :state, Set([:diffusion])) for i in 1:nvars]
        layout = VariableLayout(vars)
        SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())
    end

    @testset "FakeHotgatesModel — FD Jacobian matches analytic (via hotgates_jacobian!)" begin
        nx    = 4
        nvars = 2
        ctx   = make_ctx(nx; nvars = nvars)
        n     = nvars * nx

        fake    = FakeHotgatesModel(1.5, 0.3)
        adaptor = HotgatesTrappingAdaptor([1], [2], ["c"], ["theta"], String[], zeros(0, nx))
        op      = HotgatesReactionOperator(fake, adaptor, ConstantTemperature(300.0))

        # FakeHotgatesModel has no analytic Jacobian — uses FD fallback internally.
        # We verify the FD fallback via jacobian! matches an external FD of rhs!.
        @test !supports_jacobian(op)

        # Manually call jacobian! (the internal FD fallback is always available)
        u0 = clamp.(abs.(randn(n)) .* 0.5, 0.01, 0.99)

        f! = (du, u, ctx, t) -> begin
            fill!(du, 0.0)
            rhs!(du, op, u, ctx, t)
        end

        J_fd_ext = fd_jacobian(f!, u0, ctx)

        J_fallback = zeros(n, n)
        Flopsy.jacobian!(J_fallback, op, u0, ctx, 0.0)

        @test norm(J_fallback - J_fd_ext) / (norm(J_fd_ext) + 1e-12) < 1e-4
    end
end

@testset "Selective sparsity prototype for HotgatesReactionOperator" begin
    function make_ctx(nx; nvars=1)
        mesh = Mesh1D(1.0, nx)
        vars = [VariableInfo(Symbol("u$i"), :state, Set([:diffusion])) for i in 1:nvars]
        layout = VariableLayout(vars)
        SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())
    end

    @testset "trap_groups singleton — prototype is sparse (no cross-trap entries)" begin
        # 1 mobile + 3 independent traps (each their own group)
        nx    = 3
        nvars = 4   # 1 mobile + 3 traps
        ctx   = make_ctx(nx; nvars = nvars)
        n     = nvars * nx

        fake    = FakeHotgatesModel(1.0, 0.5)
        # 3 traps, each its own singleton group (default)
        adaptor = HotgatesTrappingAdaptor(
            [1], [2, 3, 4],
            ["c"], ["th1", "th2", "th3"],
            String[], zeros(0, nx),
        )
        # default 6-arg: trap_groups = [[1],[2],[3]] → no inter-trap coupling
        @test adaptor.trap_groups == [[1], [2], [3]]

        diff  = LinearDiffusionOperator([1.0, 0.0, 0.0, 0.0], layout -> [1], nothing)
        react = HotgatesReactionOperator(fake, adaptor, ConstantTemperature(300.0))
        ops   = (diff, react)
        model = SystemModel(ctx.layout, (;), ctx)

        proto_selective = Flopsy._build_jac_prototype(model, ops)
        nnz_selective   = nnz(proto_selective)

        # Dense diagonal blocks would be nvars² per node + diffusion off-diagonal.
        # Selective pattern excludes cross-trap entries (e.g. (th1,th2) pairs).
        # With 1 mobile + 3 independent traps and n=4 per node, full diagonal block = 4²=16.
        # Selective per-node = mobile-mobile(1) + mobile-trap(3) + trap-mobile(3) + diagonal-traps(3) = 10.
        # Plus off-diagonal diffusion (nx-1)=2 entries × 2 directions = 4 entries.
        nnz_dense = nvars^2 * nx + 2 * (nx - 1)   # fully dense diagonal + diffusion off-diag
        @test nnz_selective < nnz_dense
    end

    @testset "trap_groups coupled — tridiagonal pattern preserved in prototype" begin
        # 1 mobile + 4 traps in one group (multi-occupancy, adjacent levels coupled)
        nx    = 3
        nvars = 5   # 1 mobile + 4 trap levels
        ctx   = make_ctx(nx; nvars = nvars)
        n     = nvars * nx

        fake    = FakeHotgatesModel(1.0, 0.5)
        # One defect type with 4 levels — pass trap_groups explicitly
        adaptor = HotgatesTrappingAdaptor(
            [1], [2, 3, 4, 5],
            ["c"], ["th1", "th2", "th3", "th4"],
            String[], zeros(0, nx),
            [[1, 2, 3, 4]],   # one group of 4 → tridiagonal 4×4 trap block
        )

        diff  = LinearDiffusionOperator([1.0, 0.0, 0.0, 0.0, 0.0], layout -> [1], nothing)
        react = HotgatesReactionOperator(fake, adaptor, ConstantTemperature(300.0))
        ops   = (diff, react)
        model = SystemModel(ctx.layout, (;), ctx)

        proto = Flopsy._build_jac_prototype(model, ops)
        I_p, J_p, _ = findnz(proto)
        proto_set = Set(zip(I_p, J_p))

        # For node ix=1 (offset=0), check tridiagonal trap-trap entries are present
        for k in 1:3
            r  = 1 + k      # trap var k+1 (1-based row in full vector)
            r1 = r + 1      # next trap level
            @test (r, r)  ∈ proto_set   # diagonal
            @test (r, r1) ∈ proto_set   # super-diagonal
            @test (r1, r) ∈ proto_set   # sub-diagonal
        end

        # Corner entry (th1, th4) should NOT be in the prototype
        @test (2, 5) ∉ proto_set
        @test (5, 2) ∉ proto_set
    end
end
