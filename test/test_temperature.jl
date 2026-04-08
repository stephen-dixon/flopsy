@testset "Temperature providers" begin
    nx = 5
    mesh = Mesh1D(1.0, nx)
    layout = VariableLayout([VariableInfo(:u, :state, Set{Symbol}())])
    ctx = SystemContext(layout, nx, mesh, Dict{Symbol,Any}(), Dict{Symbol,Any}())

    @testset "ConstantTemperature" begin
        tp = ConstantTemperature(300.0)
        @test temperature_at(tp, ctx, 0.0, 1) == 300.0
        @test temperature_at(tp, ctx, 100.0, 3) == 300.0
        @test temperature_at(tp, ctx, 1e6, nx) == 300.0
    end

    @testset "LinearRampTemperature" begin
        T0 = 300.0
        rate = 10.0 / 60.0   # 10 K/min in K/s
        tp = LinearRampTemperature(T0, rate)
        @test temperature_at(tp, ctx, 0.0, 1) ≈ T0
        @test temperature_at(tp, ctx, 60.0, 1) ≈ T0 + 10.0
        @test temperature_at(tp, ctx, 600.0, 2) ≈ T0 + 100.0
        # spatially uniform — same for all nodes
        @test temperature_at(tp, ctx, 60.0, nx) ≈ temperature_at(tp, ctx, 60.0, 1)
    end

    @testset "FunctionTemperature" begin
        f = t -> 300.0 + 0.5 * t
        tp = FunctionTemperature(f)
        @test temperature_at(tp, ctx, 0.0, 1) ≈ 300.0
        @test temperature_at(tp, ctx, 100.0, 1) ≈ 350.0
        @test temperature_at(tp, ctx, 100.0, nx) ≈ 350.0  # spatially uniform
    end

    @testset "FunctionTemperature piecewise" begin
        # Simulate a loading → ramp profile
        T_load = 400.0
        T_rest = 300.0
        t_rest = 86400.0
        ramp_rate = 10.0 / 60.0
        f = t -> begin
            t < t_rest && return T_load
            T_rest + ramp_rate * (t - t_rest)
        end
        tp = FunctionTemperature(f)
        @test temperature_at(tp, ctx, 0.0, 1) ≈ T_load
        @test temperature_at(tp, ctx, t_rest - 1, 1) ≈ T_load
        @test temperature_at(tp, ctx, t_rest, 1) ≈ T_rest
        @test temperature_at(tp, ctx, t_rest + 60.0, 1) ≈ T_rest + 10.0
    end
end
