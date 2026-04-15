"""
    AbstractDiffusionCoefficients

Abstract type for objects that supply diffusion coefficients per variable,
node, and temperature.  Implement `get_D(coefficients, ivar, ix, T)`.

Used by `LinearDiffusionOperator` and `DirichletBoundaryOperator` to support
constant, Arrhenius, callable, and external-library (e.g. Palioxis) coefficients
within the same operator interface.

Concrete types:
- `ConstantDiffusion`  — fixed values independent of position and temperature
- `ArrheniusDiffusion` — D(T) = D₀ exp(−Eₐ / kBT), spatially uniform
- `CallableDiffusion`  — one callable `f(T)` per variable, spatially uniform
"""
abstract type AbstractDiffusionCoefficients end

"""
    get_D(coefficients, ivar, ix, T) -> Float64

Return the diffusion coefficient for variable index `ivar` at spatial node `ix`
and temperature `T` (Kelvin).

Dispatches on the concrete type of `coefficients`.  For `AbstractVector{<:Real}`,
the value is returned directly (constant: both `ix` and `T` are ignored).

The `ix` argument is reserved for future spatially-varying diffusion models.
All current concrete types ignore it.
"""
function get_D end

# Plain vector fallback — constant coefficients, ix and T are ignored.
get_D(coeffs::AbstractVector{<:Real}, ivar::Int, ix::Int, T) = Float64(coeffs[ivar])

# ---------------------------------------------------------------------------

"""
    ConstantDiffusion(values)

Temperature-independent, spatially-uniform diffusion coefficients.
`values[ivar]` is used directly regardless of position or temperature.

Equivalent to passing a plain `Vector{Float64}` to `LinearDiffusionOperator`, but
allows mixing with other `AbstractDiffusionCoefficients` types in operator trees.
"""
struct ConstantDiffusion{T <: Real} <: AbstractDiffusionCoefficients
    values::Vector{T}
end

get_D(c::ConstantDiffusion, ivar::Int, ix::Int, T) = Float64(c.values[ivar])

# ---------------------------------------------------------------------------

"""
    ArrheniusDiffusion(D0, Ea; kb=8.617333e-5)

Arrhenius temperature-dependent diffusion:

    D(T) = D₀[ivar] * exp(−Eₐ[ivar] / (kB * T))

- `D0` — pre-exponential factors [m²/s], one per variable
- `Ea` — activation energies [eV], one per variable
- `kb` — Boltzmann constant [eV/K] (default: 8.617333 × 10⁻⁵ eV/K)

`T` must be in Kelvin.  The coefficient is spatially uniform (same at every node).

# Example
```julia
# Single mobile species: D₀ = 1e-7 m²/s, Eₐ = 0.2 eV
D_mobile = ArrheniusDiffusion([1e-7], [0.2])
diffusion = LinearDiffusionOperator(D_mobile, selector, nothing, temperature)
```
"""
struct ArrheniusDiffusion{T <: Real} <: AbstractDiffusionCoefficients
    D0::Vector{T}
    Ea::Vector{T}
    kb::Float64
end

function ArrheniusDiffusion(D0::AbstractVector{T}, Ea::AbstractVector{T};
        kb::Float64 = 8.617333e-5) where {T <: Real}
    length(D0) == length(Ea) || throw(ArgumentError("D0 and Ea must have the same length"))
    ArrheniusDiffusion(collect(T, D0), collect(T, Ea), kb)
end

function get_D(c::ArrheniusDiffusion, ivar::Int, ix::Int, T::Real)
    c.D0[ivar] * exp(-c.Ea[ivar] / (c.kb * T))
end

# ---------------------------------------------------------------------------

"""
    CallableDiffusion(funcs)

Diffusion coefficients given by one callable per variable.  Each element of `funcs`
must be a callable `f(T::Real) -> Float64`.

Useful for arbitrary analytic expressions or lookup-table interpolations.
The coefficient is spatially uniform (the callable receives only `T`, not `ix`).

# Example
```julia
# Two species: constant D for mobile, zero for trapped
D_mobile = T -> 1e-7 * exp(-0.2 / (8.617e-5 * T))
D_trap   = T -> 0.0
diffusion_coeffs = CallableDiffusion([D_mobile, D_trap])
```
"""
struct CallableDiffusion{F} <: AbstractDiffusionCoefficients
    funcs::Vector{F}
end

get_D(c::CallableDiffusion, ivar::Int, ix::Int, T::Real) = Float64(c.funcs[ivar](T))
