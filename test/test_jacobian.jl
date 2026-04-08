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

    @testset "DirichletBoundaryOperator — both BCs" begin
        nx  = 6
        ctx = make_ctx(nx)
        D   = 1.0
        op  = DirichletBoundaryOperator(selector_all, [D], nothing;
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
