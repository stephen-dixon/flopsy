"""
    AbstractDiffusionCoefficients

Abstract type for objects that supply diffusion coefficients per variable and
temperature.  Implement `get_D(coefficients, ivar, T)`.

Used by `LinearDiffusionOperator` and `DirichletBoundaryOperator` to support
constant, Arrhenius, callable, and external-library (e.g. Palioxis) coefficients
within the same operator interface.

Concrete types:
- `ConstantDiffusion`  — fixed values independent of temperature
- `ArrheniusDiffusion` — D(T) = D₀ exp(−Eₐ / kBT)
- `CallableDiffusion`  — one callable `f(T)` per variable
"""
abstract type AbstractDiffusionCoefficients end


"""
    get_D(coefficients, ivar, T) -> Float64

Return the diffusion coefficient for variable index `ivar` at temperature `T` (Kelvin).

Dispatches on the concrete type of `coefficients`.  For `AbstractVector{<:Real}`,
the value is returned directly (constant, T is ignored).
"""
function get_D end

# Plain vector fallback — constant coefficients, T is ignored.
get_D(coeffs::AbstractVector{<:Real}, ivar::Int, T) = Float64(coeffs[ivar])


# ---------------------------------------------------------------------------

"""
    ConstantDiffusion(values)

Temperature-independent diffusion coefficients.  `values[ivar]` is used directly.

Equivalent to passing a plain `Vector{Float64}` to `LinearDiffusionOperator`, but
allows mixing with other `AbstractDiffusionCoefficients` types in operator trees.
"""
struct ConstantDiffusion{T<:Real} <: AbstractDiffusionCoefficients
    values::Vector{T}
end

get_D(c::ConstantDiffusion, ivar::Int, T) = Float64(c.values[ivar])


# ---------------------------------------------------------------------------

"""
    ArrheniusDiffusion(D0, Ea; kb=8.617333e-5)

Arrhenius temperature-dependent diffusion:

    D(T) = D₀[ivar] * exp(−Eₐ[ivar] / (kB * T))

- `D0` — pre-exponential factors [m²/s], one per variable
- `Ea` — activation energies [eV], one per variable
- `kb` — Boltzmann constant [eV/K] (default: 8.617333 × 10⁻⁵ eV/K)

`T` must be in Kelvin.

# Example
```julia
# Single mobile species: D₀ = 1e-7 m²/s, Eₐ = 0.2 eV
D_mobile = ArrheniusDiffusion([1e-7], [0.2])
diffusion = LinearDiffusionOperator(D_mobile, selector, nothing, temperature)
```
"""
struct ArrheniusDiffusion{T<:Real} <: AbstractDiffusionCoefficients
    D0::Vector{T}
    Ea::Vector{T}
    kb::Float64
end

function ArrheniusDiffusion(D0::AbstractVector{T}, Ea::AbstractVector{T};
                             kb::Float64=8.617333e-5) where {T<:Real}
    length(D0) == length(Ea) || throw(ArgumentError("D0 and Ea must have the same length"))
    ArrheniusDiffusion(collect(T, D0), collect(T, Ea), kb)
end

get_D(c::ArrheniusDiffusion, ivar::Int, T::Real) = c.D0[ivar] * exp(-c.Ea[ivar] / (c.kb * T))


# ---------------------------------------------------------------------------

"""
    CallableDiffusion(funcs)

Diffusion coefficients given by one callable per variable.  Each element of `funcs`
must be a callable `f(T::Real) -> Float64`.

Useful for arbitrary analytic expressions or lookup-table interpolations.

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

get_D(c::CallableDiffusion, ivar::Int, T::Real) = Float64(c.funcs[ivar](T))
