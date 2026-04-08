"""
    AbstractTemperatureProvider

Abstract type for objects that supply a temperature value at a given spatial
node and time.  Implement `temperature_at(provider, ctx, t, ix)`.
"""
abstract type AbstractTemperatureProvider end


"""
    ConstantTemperature(value)

Returns the same temperature at every node and time.
"""
struct ConstantTemperature{T} <: AbstractTemperatureProvider
    value::T
end

temperature_at(tp::ConstantTemperature, ctx, t, ix) = tp.value


"""
    LinearRampTemperature(T0, ramp_rate)

Temperature that increases linearly in time:

    T(t) = T0 + ramp_rate * t

`T0` is the initial temperature in Kelvin; `ramp_rate` is in K/s.
Typical TDS ramp rates are 0.5–5 K/s.
"""
struct LinearRampTemperature{T} <: AbstractTemperatureProvider
    T0::T
    ramp_rate::T
end

temperature_at(tp::LinearRampTemperature, ctx, t, ix) = tp.T0 + tp.ramp_rate * t


"""
    FunctionTemperature(f)

Temperature given by an arbitrary callable `f(t)` returning Kelvin.
The same value is used at every spatial node (uniform temperature field).

Example:
```julia
# Ramp from 300 K to 1000 K over 700 s, then hold
T_func = FunctionTemperature(t -> min(300.0 + t, 1000.0))
```
"""
struct FunctionTemperature{F} <: AbstractTemperatureProvider
    f::F
end

temperature_at(tp::FunctionTemperature, ctx, t, ix) = tp.f(t)
