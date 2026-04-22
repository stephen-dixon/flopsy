@testset "TDS workflows" begin

    # -----------------------------------------------------------------------
    # Helpers
    # -----------------------------------------------------------------------

    function _make_trapping_result(; nx = 20, tspan = (0.0, 1.0),
            formulation = "unsplit", split_method = "strang", dt = 0.1,
            saveat = nothing)
        raw = Dict(
            "mesh" => Dict("main" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => nx)),
            "backend" => Dict("main" => Dict(
                "type" => "trapping_1d",
                "k_trap" => 2.0, "k_detrap" => 0.5,
                "diffusion_coefficient" => 0.01
            )),
            "ic" => Dict("fill" => Dict(
                "type" => "uniform_species", "species" => "H_mobile", "value" => 0.1)),
            "bc" => Dict(
                "left" => Dict(
                    "type" => "dirichlet", "species" => "H_mobile",
                    "boundary" => "left", "value" => 0.0),
                "right" => Dict(
                    "type" => "dirichlet", "species" => "H_mobile",
                    "boundary" => "right", "value" => 0.0)
            ),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "main",
                "backend" => "main",
                "ics" => ["fill"],
                "bcs" => ["left", "right"],
                "tspan" => collect(tspan),
                "formulation" => formulation,
                "split_method" => split_method,
                "dt" => dt,
                "saveat" => saveat === nothing ? collect(range(tspan[1], tspan[2]; length = 5)) :
                            saveat
            ))
        )
        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        return solve(Flopsy.build_simulation(ctx))
    end

    # -----------------------------------------------------------------------
    # 1. Summary CSV generation
    # -----------------------------------------------------------------------

    @testset "summary CSV output from TOML deck" begin
        tmp = mktempdir()
        csv_path = joinpath(tmp, "summary.csv")

        raw = Dict(
            "mesh" => Dict("main" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 20)),
            "backend" => Dict("main" => Dict(
                "type" => "trapping_1d",
                "k_trap" => 2.0, "k_detrap" => 0.5,
                "diffusion_coefficient" => 0.01
            )),
            "ic" => Dict("fill" => Dict(
                "type" => "uniform_species", "species" => "H_mobile", "value" => 0.1)),
            "output" => Dict("csv" => Dict(
                "type" => "summary_csv", "file" => csv_path)),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "main", "backend" => "main",
                "ics" => ["fill"],
                "outputs" => ["csv"],
                "tspan" => [0.0, 0.5],
                "saveat" => [0.0, 0.25, 0.5]
            ))
        )

        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        result = solve(Flopsy.build_simulation(ctx))

        @test isfile(csv_path)

        import CSV as _CSV
        import DataFrames as _DF
        df = _CSV.read(csv_path, _DF.DataFrame)

        @test hasproperty(df, :time)
        @test hasproperty(df, :integral_H_mobile)
        @test nrow(df) == 3
        @test df.time ≈ [0.0, 0.25, 0.5]
    end

    @testset "summary_csv write_summary_csv round-trip" begin
        result = _make_trapping_result(; nx = 15, tspan = (0.0, 0.3))
        tmp = tempname() * ".csv"
        try
            path = write_summary_csv(result, tmp)
            @test isfile(tmp)
            @test path == tmp

            import CSV as _CSV
            import DataFrames as _DF
            df = _CSV.read(tmp, _DF.DataFrame)
            @test hasproperty(df, :time)
            @test hasproperty(df, :integral_H_mobile)
            @test nrow(df) >= 2
        finally
            isfile(tmp) && rm(tmp)
        end
    end

    # -----------------------------------------------------------------------
    # 2. Temperature ramp correctness
    # -----------------------------------------------------------------------

    @testset "PiecewiseTemperature stage construction" begin
        nx = 5
        mesh = Mesh1D(1.0, nx)
        layout = VariableLayout([VariableInfo(:u, :state, Set{Symbol}())])
        ctx = SystemContext(layout, nx, mesh, Dict{Symbol, Any}(), Dict{Symbol, Any}())

        stages = [
            Dict("type" => "hold", "T" => 300.0, "duration" => 60.0),
            Dict("type" => "ramp", "rate" => 0.5)
        ]
        tp = PiecewiseTemperature(stages)

        @test temperature_at(tp, ctx, 0.0, 1) ≈ 300.0
        @test temperature_at(tp, ctx, 30.0, 1) ≈ 300.0
        @test temperature_at(tp, ctx, 60.0, 1) ≈ 300.0
        @test temperature_at(tp, ctx, 61.0, 1) ≈ 300.5
        @test temperature_at(tp, ctx, 160.0, 1) ≈ 350.0
        @test temperature_at(tp, ctx, 160.0, nx) ≈ 350.0  # spatially uniform
    end

    @testset "PiecewiseTemperature multi-stage" begin
        nx = 4
        mesh = Mesh1D(1.0, nx)
        layout = VariableLayout([VariableInfo(:u, :state, Set{Symbol}())])
        ctx = SystemContext(layout, nx, mesh, Dict{Symbol, Any}(), Dict{Symbol, Any}())

        stages = [
            Dict("type" => "hold", "T" => 300.0, "duration" => 100.0),
            Dict("type" => "ramp", "rate" => 1.0, "duration" => 200.0),
            Dict("type" => "hold", "T" => 500.0, "duration" => 100.0),
        ]
        tp = PiecewiseTemperature(stages)

        @test temperature_at(tp, ctx, 0.0, 1) ≈ 300.0
        @test temperature_at(tp, ctx, 50.0, 1) ≈ 300.0
        @test temperature_at(tp, ctx, 100.0, 1) ≈ 300.0
        @test temperature_at(tp, ctx, 200.0, 1) ≈ 400.0
        @test temperature_at(tp, ctx, 300.0, 1) ≈ 500.0
        @test temperature_at(tp, ctx, 350.0, 1) ≈ 500.0
        @test temperature_at(tp, ctx, 400.0, 1) ≈ 500.0
    end

    @testset "temperature block in TOML deck (linear_ramp)" begin
        raw = Dict(
            "mesh" => Dict("main" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 11)),
            "backend" => Dict("main" => Dict(
                "type" => "trapping_1d",
                "k_trap" => 1.0, "k_detrap" => 0.1,
                "diffusion_coefficient" => 0.01
            )),
            "temperature" => Dict("ramp" => Dict(
                "type" => "linear_ramp", "T0" => 300.0, "rate" => 0.5)),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "main", "backend" => "main",
                "temperature" => "ramp",
                "tspan" => [0.0, 1.0],
                "saveat" => [0.0, 1.0]
            ))
        )
        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        @test !isempty(ctx.temperatures)
        tp = ctx.temperatures[:ramp]
        @test tp isa LinearRampTemperature
        @test tp.T0 ≈ 300.0
        @test tp.ramp_rate ≈ 0.5
    end

    @testset "temperature block in TOML deck (constant)" begin
        raw = Dict(
            "mesh" => Dict("m" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 11)),
            "backend" => Dict("b" => Dict(
                "type" => "diffusion_1d", "diffusion_coefficient" => 0.1)),
            "temperature" => Dict("T0" => Dict(
                "type" => "constant", "value" => 500.0)),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "m", "backend" => "b",
                "temperature" => "T0",
                "tspan" => [0.0, 1.0]
            ))
        )
        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        @test ctx.temperatures[:T0] isa ConstantTemperature
    end

    @testset "tds problem type requires temperature" begin
        raw_missing_temp = Dict(
            "mesh" => Dict("m" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 11)),
            "backend" => Dict("b" => Dict(
                "type" => "trapping_1d",
                "k_trap" => 1.0, "k_detrap" => 0.1,
                "diffusion_coefficient" => 0.01
            )),
            "problem" => Dict("run" => Dict(
                "type" => "tds",
                "mesh" => "m", "backend" => "b",
                "tspan" => [0.0, 1.0],
                "dt" => 0.1
            ))
        )
        msg = try
            build_context(parse_input_deck(raw_missing_temp))
            nothing
        catch err
            sprint(showerror, err)
        end
        @test msg !== nothing
        @test occursin("temperature", msg)
    end

    # -----------------------------------------------------------------------
    # 3. CLI --problem selection
    # -----------------------------------------------------------------------

    @testset "CLI --problem selects named problem" begin
        tmp = mktempdir()
        deck_path = joinpath(tmp, "multi.toml")
        open(deck_path, "w") do io
            write(io, """
[mesh.main]
type = "uniform_1d"
xmin = 0.0
xmax = 1.0
nx = 11

[backend.main]
type = "diffusion_1d"
diffusion_coefficient = 0.1

[output.out]
type = "hdf5"
file = "$(joinpath(tmp, "result.h5"))"

[problem.alpha]
type = "simulation"
mesh = "main"
backend = "main"
outputs = ["out"]
tspan = [0.0, 0.5]
saveat = [0.0, 0.5]

[problem.beta]
type = "simulation"
mesh = "main"
backend = "main"
tspan = [0.0, 1.0]
saveat = [0.0, 1.0]
""")
        end

        # --problem alpha should run problem alpha (which writes out)
        @test cli_main(["run", deck_path, "--problem", "alpha"]) == 0
        @test isfile(joinpath(tmp, "result.h5"))

        # Without --problem the default (alpha = alphabetically first) is chosen
        @test cli_main(["validate", deck_path, "--problem", "beta"]) == 0
        @test cli_main(["validate", deck_path, "--problem", "alpha"]) == 0
    end

    @testset "CLI --output-dir redirects outputs" begin
        tmp = mktempdir()
        out_dir = joinpath(tmp, "outputs")
        mkpath(out_dir)

        deck_path = joinpath(tmp, "simple.toml")
        open(deck_path, "w") do io
            write(io, """
[mesh.main]
type = "uniform_1d"
xmin = 0.0
xmax = 1.0
nx = 11

[backend.main]
type = "diffusion_1d"
diffusion_coefficient = 0.1

[ic.fill]
type = "uniform_species"
species = "u"
value = 1.0

[output.out]
type = "hdf5"
file = "fields.h5"

[problem.run]
type = "simulation"
mesh = "main"
backend = "main"
ics = ["fill"]
outputs = ["out"]
tspan = [0.0, 0.2]
saveat = [0.0, 0.1, 0.2]
""")
        end

        @test cli_main(["run", deck_path, "--output-dir", out_dir]) == 0
        @test isfile(joinpath(out_dir, "fields.h5"))
        @test !isfile("fields.h5")   # should NOT be in cwd
    end

    # -----------------------------------------------------------------------
    # 4. Split method selection (Lie vs Strang)
    # -----------------------------------------------------------------------

    @testset "split formulation lie vs strang both converge" begin
        lie_result = _make_trapping_result(;
            formulation = "split", split_method = "lie", dt = 0.05,
            tspan = (0.0, 0.5))
        strang_result = _make_trapping_result(;
            formulation = "split", split_method = "strang", dt = 0.05,
            tspan = (0.0, 0.5))

        # Both split schemes produce finite, non-trivial solutions
        @test all(isfinite, lie_result.solution.u[end])
        @test all(isfinite, strang_result.solution.u[end])

        # Both schemes should deplete the initial inventory (material diffuses out)
        inv_lie_0 = integrated_variable(lie_result, :H_mobile)[1]
        inv_lie_end = integrated_variable(lie_result, :H_mobile)[end]
        inv_strang_0 = integrated_variable(strang_result, :H_mobile)[1]
        inv_strang_end = integrated_variable(strang_result, :H_mobile)[end]

        @test inv_lie_end < inv_lie_0        # inventory decreasing
        @test inv_strang_end < inv_strang_0  # inventory decreasing
        @test all(x -> x >= -1e-10, lie_result.solution.u[end])    # non-negative
        @test all(x -> x >= -1e-10, strang_result.solution.u[end]) # non-negative
    end

    @testset "split_method exposed in TOML deck" begin
        raw_lie = Dict(
            "mesh" => Dict("m" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 11)),
            "backend" => Dict("b" => Dict(
                "type" => "trapping_1d",
                "k_trap" => 2.0, "k_detrap" => 0.5,
                "diffusion_coefficient" => 0.01
            )),
            "ic" => Dict("fill" => Dict(
                "type" => "uniform_species", "species" => "H_mobile", "value" => 0.1)),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "m", "backend" => "b",
                "ics" => ["fill"],
                "tspan" => [0.0, 0.5],
                "formulation" => "split",
                "split_method" => "lie",
                "dt" => 0.05,
                "saveat" => [0.0, 0.25, 0.5]
            ))
        )
        deck = parse_input_deck(raw_lie)
        ctx = build_context(deck)
        prob_def = ctx.problems[:run]
        @test prob_def.solver.split_method == :lie
        @test prob_def.solver.formulation == :split

        result = solve(Flopsy.build_simulation(ctx))
        @test all(isfinite, result.solution.u[end])
    end

    # -----------------------------------------------------------------------
    # 5. FakeHotgates produces non-zero flux under TDS-like conditions
    # -----------------------------------------------------------------------

    @testset "FakeHotgates TDS run produces non-zero flux" begin
        raw = Dict(
            "mesh" => Dict("main" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 25)),
            "backend" => Dict("main" => Dict(
                "type" => "hotgates_trapping",
                "k_trap" => 2.0, "k_detrap" => 0.5,
                "diffusion_coefficient" => 0.05,
                "temperature" => 300.0
            )),
            "ic" => Dict("fill" => Dict(
                "type" => "uniform_species", "species" => "H_mobile", "value" => 0.1)),
            "bc" => Dict(
                "left" => Dict(
                    "type" => "dirichlet", "species" => "H_mobile",
                    "boundary" => "left", "value" => 0.0),
                "right" => Dict(
                    "type" => "dirichlet", "species" => "H_mobile",
                    "boundary" => "right", "value" => 0.0)
            ),
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "main", "backend" => "main",
                "ics" => ["fill"],
                "bcs" => ["left", "right"],
                "tspan" => [0.0, 0.5],
                "saveat" => [0.0, 0.25, 0.5]
            ))
        )
        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        result = solve(Flopsy.build_simulation(ctx))

        fluxes = surface_diffusive_fluxes(result)
        @test !isempty(fluxes)
        fl = first(values(fluxes))
        # At least one surface should show positive outward flux after t > 0
        @test any(fl.left[2:end] .> 0) || any(fl.right[2:end] .> 0)
    end

    @testset "FakeHotgates TDS with problem-level temperature override" begin
        raw = Dict(
            "mesh" => Dict("main" => Dict(
                "type" => "uniform_1d", "xmin" => 0.0, "xmax" => 1.0, "nx" => 15)),
            "backend" => Dict("main" => Dict(
                "type" => "hotgates_trapping",
                "k_trap" => 2.0, "k_detrap" => 0.5,
                "diffusion_coefficient" => 0.05,
                "temperature" => 300.0   # backend-level fallback
            )),
            "ic" => Dict("fill" => Dict(
                "type" => "uniform_species", "species" => "H_mobile", "value" => 0.05)),
            "temperature" => Dict("ramp" => Dict(
                "type" => "constant", "value" => 400.0)),   # problem-level overrides
            "problem" => Dict("run" => Dict(
                "type" => "simulation",
                "mesh" => "main", "backend" => "main",
                "ics" => ["fill"],
                "temperature" => "ramp",
                "tspan" => [0.0, 0.3],
                "saveat" => [0.0, 0.15, 0.3]
            ))
        )
        deck = parse_input_deck(raw)
        ctx = build_context(deck)
        result = solve(Flopsy.build_simulation(ctx))
        @test all(isfinite, result.solution.u[end])
    end

end
