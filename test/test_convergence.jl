@testset "Spatial convergence — diffusion O(dx²)" begin
    #=
    Convergence study for LinearDiffusionOperator applied to a known smooth function.

    We evaluate the discrete diffusion operator on u(x) = sin(π x) and compare
    the output against the exact Laplacian D * u_xx(x) = -D π² sin(π x).

    Interior nodes use the standard centered-difference stencil, which is second-order
    accurate.  Boundary nodes use a one-sided stencil (first-order), so they are
    excluded from the L2 norm — the convergence test targets the interior accuracy.

    Expected: ||LDO(u) - D u_xx||_{L2,interior} = O(dx²), so the observed order
    between successive grids should be ≥ 1.8.
    =#

    D      = 1.0
    u_fn(x)    = sin(π * x)
    uxx_fn(x)  = -π^2 * sin(π * x)

    function run_grid(nx)
        mesh   = Mesh1D(1.0, nx)
        x      = mesh.x
        dx     = mesh.dx
        vars   = [VariableInfo(:u, :state, Set([:diffusion]))]
        layout = VariableLayout(vars)
        ctx    = SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())

        op = LinearDiffusionOperator([D], layout -> [1], nothing)

        u  = u_fn.(x)
        du = zeros(nx)
        rhs!(du, op, u, ctx, 0.0)   # du ≈ D * u_xx at each node

        # Restrict to interior nodes — boundary uses first-order one-sided stencil
        interior   = 2:(nx - 1)
        du_exact   = D .* uxx_fn.(x[interior])
        l2_err     = sqrt(dx * sum((du[interior] .- du_exact).^2))
        return dx, l2_err
    end

    # Grids: 16, 32, 64, 128, 256 nodes
    grids  = [16, 32, 64, 128, 256]
    errors = [run_grid(nx) for nx in grids]
    dxs    = [e[1] for e in errors]
    l2s    = [e[2] for e in errors]

    # Estimate convergence order between consecutive grid pairs
    orders = [log(l2s[i] / l2s[i+1]) / log(dxs[i] / dxs[i+1]) for i in 1:length(grids)-1]

    # Interior stencil is O(dx²): all observed orders should be ≥ 1.8
    @test all(orders .≥ 1.8)

    # Absolute accuracy on the finest grid
    @test l2s[end] < 1e-4
end

@testset "Spatial convergence — multi-variable diffusion (selective variable)" begin
    #=
    Verify that when only one of two variables is marked for diffusion, the LDO
    applies to that variable and leaves the other unchanged.  Also confirms
    O(dx²) convergence for the diffusing variable.
    =#

    D = 2.0
    u_fn(x)   = sin(2π * x)
    uxx_fn(x) = -4π^2 * sin(2π * x)

    function run_grid(nx)
        mesh   = Mesh1D(1.0, nx)
        x      = mesh.x
        dx     = mesh.dx
        # Two variables; only var 1 diffuses
        vars   = [VariableInfo(:u1, :state, Set([:diffusion])),
                  VariableInfo(:u2, :state, Set([:nodiffusion]))]
        layout = VariableLayout(vars)
        ctx    = SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())

        op = LinearDiffusionOperator([D, 0.0], layout -> [1], nothing)

        # State: [u1[1], u2[1], u1[2], u2[2], ...] (nvars=2 layout)
        u  = zeros(2 * nx)
        dU = state_view(u, layout, nx)
        for ix in 1:nx
            dU[1, ix] = u_fn(x[ix])
            dU[2, ix] = 1.0  # constant — should be unchanged
        end

        du = zeros(2 * nx)
        rhs!(du, op, u, ctx, 0.0)

        dU_out = state_view(du, layout, nx)

        # Variable 2 must be untouched
        @test all(dU_out[2, :] .== 0.0)

        # Variable 1 interior: O(dx²)
        interior  = 2:(nx - 1)
        du1_exact = D .* uxx_fn.(x[interior])
        l2_err    = sqrt(dx * sum((dU_out[1, interior] .- du1_exact).^2))
        return dx, l2_err
    end

    grids  = [16, 32, 64, 128]
    errors = [run_grid(nx) for nx in grids]
    dxs    = [e[1] for e in errors]
    l2s    = [e[2] for e in errors]

    orders = [log(l2s[i] / l2s[i+1]) / log(dxs[i] / dxs[i+1]) for i in 1:length(grids)-1]
    @test all(orders .≥ 1.8)
end
