using Test
using Flopsy
using LinearAlgebra
using SparseArrays
using ADTypes: AutoFiniteDiff
using DataFrames: nrow
using SciMLBase: ReturnCode
using OrdinaryDiffEq: Rodas5, Rodas5P, KenCarp4

@testset "Flopsy" begin
    include("test_temperature.jl")
    include("test_diffusion_coefficients.jl")
    include("test_operators.jl")
    include("test_jacobian.jl")
    include("test_convergence.jl")
    include("test_output.jl")
    include("test_trapping.jl")
    include("test_formulations.jl")
    include("test_config.jl")
    include("test_registry_cli.jl")
    include("test_tds.jl")
end
