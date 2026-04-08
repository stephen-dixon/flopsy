using Test
using Flopsy
using LinearAlgebra
using SparseArrays
using ADTypes: AutoFiniteDiff
using DataFrames: nrow
using SciMLBase: ReturnCode
using OrdinaryDiffEq: Rodas5

@testset "Flopsy" begin
    include("test_temperature.jl")
    include("test_diffusion_coefficients.jl")
    include("test_operators.jl")
    include("test_jacobian.jl")
    include("test_output.jl")
    include("test_trapping.jl")
end
