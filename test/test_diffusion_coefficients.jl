@testset "Diffusion coefficients" begin

    @testset "Plain vector (constant)" begin
        coeffs = [1.0, 0.0, 2.5]
        @test get_D(coeffs, 1, 1, 300.0) == 1.0
        @test get_D(coeffs, 2, 3, 1000.0) == 0.0
        @test get_D(coeffs, 3, 1, NaN) == 2.5   # T ignored
    end

    @testset "ConstantDiffusion" begin
        cd = ConstantDiffusion([1.5, 0.0])
        @test get_D(cd, 1, 1, 300.0) == 1.5
        @test get_D(cd, 2, 5, 999.0) == 0.0
        @test get_D(cd, 1, 1, NaN) == 1.5       # T ignored
        @test get_D(cd, 1, 99, 300.0) == 1.5    # ix ignored
    end

    @testset "ArrheniusDiffusion" begin
        kb = 8.617333e-5
        D0 = [1e-7]
        Ea = [0.2]
        ad = ArrheniusDiffusion(D0, Ea)

        T = 600.0
        expected = D0[1] * exp(-Ea[1] / (kb * T))
        @test get_D(ad, 1, 1, T) ≈ expected rtol=1e-10

        # Higher T → larger D
        @test get_D(ad, 1, 1, 1200.0) > get_D(ad, 1, 1, 600.0)

        # ix is ignored
        @test get_D(ad, 1, 42, T) ≈ expected

        # Length mismatch raises
        @test_throws ArgumentError ArrheniusDiffusion([1e-7, 2e-7], [0.2])

        # Multiple variables
        ad2 = ArrheniusDiffusion([1e-7, 2e-8], [0.2, 0.4])
        @test get_D(ad2, 1, 1, T) ≈ 1e-7 * exp(-0.2 / (kb * T))
        @test get_D(ad2, 2, 1, T) ≈ 2e-8 * exp(-0.4 / (kb * T))
    end

    @testset "CallableDiffusion" begin
        f1 = T -> 1.5 * T
        f2 = T -> 0.0
        cd = CallableDiffusion([f1, f2])

        @test get_D(cd, 1, 1, 10.0) ≈ 15.0
        @test get_D(cd, 2, 1, 999.0) == 0.0
        # ix is ignored
        @test get_D(cd, 1, 7, 10.0) ≈ 15.0
    end
end
