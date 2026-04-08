@testset "Output" begin

    # Minimal fake solution for testing post-processing without running the solver.
    struct FakeSol
        t::Vector{Float64}
        u::Vector{Vector{Float64}}
    end

    function make_diffusion_result(; nx=8, D=1.0)
        mesh   = Mesh1D(1.0, nx)
        vars   = [VariableInfo(:c, :mobile, Set([:diffusion]))]
        layout = VariableLayout(vars)
        sel    = layout -> variables_with_tag(layout, :diffusion)
        diffop = LinearDiffusionOperator([D], sel, nothing)
        model  = build_rd_model(layout=layout, mesh=mesh, diffusion=diffop)

        # Two time steps: uniform interior, zero at boundaries
        u0 = zeros(nx); u0[2:nx-1] .= 1.0
        u1 = u0 .* 0.9
        sol = FakeSol([0.0, 0.1], [u0, u1])
        return wrap_result(model, sol, Dict{String,Any}())
    end

    @testset "surface_diffusive_fluxes sign convention" begin
        result = make_diffusion_result(nx=8, D=1.5)
        fluxes = surface_diffusive_fluxes(result)

        @test haskey(fluxes, :c)
        fl = fluxes[:c]

        # With interior=1, boundary=0:
        #   left_flux  = D*(U[2]-U[1])/dx > 0 (material leaving left surface)
        #   right_flux = D*(U[nx-1]-U[nx])/dx > 0 (material leaving right surface)
        @test fl.left[1]  > 0.0
        @test fl.right[1] > 0.0
        @test length(fl.left) == 2
        @test length(fl.right) == 2
    end

    @testset "surface_diffusive_fluxes — known value" begin
        nx  = 5
        D   = 2.0
        dx  = 1.0 / (nx - 1)

        mesh   = Mesh1D(1.0, nx)
        vars   = [VariableInfo(:c, :mobile, Set([:diffusion]))]
        layout = VariableLayout(vars)
        sel    = layout -> variables_with_tag(layout, :diffusion)
        diffop = LinearDiffusionOperator([D], sel, nothing)
        model  = build_rd_model(layout=layout, mesh=mesh, diffusion=diffop)

        # Linear profile: U[ix] = ix - 1 (so U[1]=0, U[nx]=nx-1)
        u0  = Float64.(0:nx-1)
        sol = FakeSol([0.0], [u0])
        result = wrap_result(model, sol, Dict{String,Any}())

        fluxes = surface_diffusive_fluxes(result)
        fl     = fluxes[:c]

        expected_left  = D * (u0[2] - u0[1]) / dx
        expected_right = D * (u0[nx-1] - u0[nx]) / dx

        @test fl.left[1]  ≈ expected_left  atol=1e-12
        @test fl.right[1] ≈ expected_right atol=1e-12
    end

    @testset "surface_diffusive_fluxes — no diffusion operator" begin
        nx     = 5
        mesh   = Mesh1D(1.0, nx)
        vars   = [VariableInfo(:c, :mobile, Set{Symbol}())]
        layout = VariableLayout(vars)
        model  = build_rd_model(layout=layout, mesh=mesh)

        u0  = zeros(nx)
        sol = FakeSol([0.0], [u0])
        result = wrap_result(model, sol, Dict{String,Any}())

        fluxes = surface_diffusive_fluxes(result)
        @test isempty(fluxes)
    end

    @testset "surface_diffusive_fluxes — T-dependent D" begin
        nx   = 6
        D0   = [1e-3]
        Ea   = [0.1]
        kb   = 8.617333e-5
        T    = 500.0

        mesh    = Mesh1D(1.0, nx)
        vars    = [VariableInfo(:c, :mobile, Set([:diffusion]))]
        layout  = VariableLayout(vars)
        sel     = layout -> variables_with_tag(layout, :diffusion)
        tp      = ConstantTemperature(T)
        diffop  = LinearDiffusionOperator(ArrheniusDiffusion(D0, Ea), sel, nothing, tp)
        model   = build_rd_model(layout=layout, mesh=mesh, diffusion=diffop)

        dx   = 1.0 / (nx - 1)
        u0   = zeros(nx); u0[2:nx-1] .= 1.0
        sol  = FakeSol([0.0], [u0])
        result = wrap_result(model, sol, Dict{String,Any}())

        fluxes = surface_diffusive_fluxes(result)
        D_expected = D0[1] * exp(-Ea[1] / (kb * T))
        @test fluxes[:c].left[1]  ≈ D_expected * (u0[2] - u0[1]) / dx atol=1e-12
    end

    @testset "HDF5 round-trip via load_ic_from_hdf5" begin
        nx   = 8
        mesh = Mesh1D(1.0, nx)
        vars = [VariableInfo(:c, :mobile, Set([:diffusion])),
                VariableInfo(:theta, :trap, Set{Symbol}())]
        layout = VariableLayout(vars)
        sel    = layout -> variables_with_tag(layout, :diffusion)
        react  = SimpleTrappingReactionOperator(2.0, 0.5, 1, 2)
        diffop = LinearDiffusionOperator([1.0, 0.0], sel, nothing)
        model  = build_rd_model(layout=layout, mesh=mesh, reaction=react, diffusion=diffop)

        # Run a short simulation to get a non-trivial solution
        u0 = zeros(2 * nx)
        U0 = state_view(u0, layout, nx)
        U0[1, :] .= 0.3   # initial mobile
        U0[2, :] .= 0.1   # initial trap occupancy

        cfg = SolverConfig(
            formulation = UnsplitFormulation(),
            algorithm   = Rodas5(autodiff=AutoFiniteDiff()),
            saveat      = [0.0, 0.05, 0.1],
        )

        sol    = solve_problem(model, u0, (0.0, 0.1), cfg)
        result = wrap_result(model, sol, Dict{String,Any}())

        tmp = tempname() * ".h5"
        try
            write_field_output_hdf5(result, tmp)

            # Load last time step as IC
            u_loaded = load_ic_from_hdf5(tmp, model)

            # Should match the final saved state
            u_final = result.solution.u[end]
            @test u_loaded ≈ u_final rtol=1e-10

            # Load first time step
            u_first = load_ic_from_hdf5(tmp, model; time_index=1)
            @test u_first ≈ result.solution.u[1] rtol=1e-10
        finally
            isfile(tmp) && rm(tmp)
        end
    end

    @testset "build_summary_dataframe columns" begin
        result = make_diffusion_result(nx=6, D=1.0)
        df = build_summary_dataframe(result)

        @test hasproperty(df, :time)
        @test hasproperty(df, :integral_c)
        @test hasproperty(df, :left_flux_c)
        @test hasproperty(df, :right_flux_c)
        @test nrow(df) == 2
    end

    @testset "library_versions contains flopsy_version" begin
        vs = Flopsy.library_versions()
        @test haskey(vs, "flopsy_version")
        @test !isempty(vs["flopsy_version"])
        # Should be a valid version string (not the old hardcoded "0.1.0-dev")
        @test vs["flopsy_version"] != "0.1.0-dev"
    end
end
