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

"""
    temperature_at(provider, ctx, t, ix) -> Real

Return the temperature in Kelvin at node `ix` and time `t`.
"""
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

"""
    PiecewiseTemperature(breakpoints, T_at_breakpoints, rates)

Piecewise-linear temperature profile built from explicit segment data.
Each segment `i` spans `[breakpoints[i], breakpoints[i+1])` with
`T(t) = T_at_breakpoints[i] + rates[i] * (t - breakpoints[i])`.
The last segment extends to infinity.

Use the stage-based constructor for input-deck-driven construction.
"""
struct PiecewiseTemperature <: AbstractTemperatureProvider
    breakpoints::Vector{Float64}
    T_at_breakpoints::Vector{Float64}
    rates::Vector{Float64}
end

function temperature_at(tp::PiecewiseTemperature, ctx, t, ix)
    n = length(tp.breakpoints)
    n == 0 && return 0.0
    for i in n:-1:1
        t >= tp.breakpoints[i] && return tp.T_at_breakpoints[i] + tp.rates[i] * (t - tp.breakpoints[i])
    end
    return tp.T_at_breakpoints[1]
end

"""
    PiecewiseTemperature(stages)

Build a `PiecewiseTemperature` from a vector of stage dicts (or NamedTuples).

Each stage must have a `"type"` key (`"hold"` or `"ramp"`).

- `"hold"`: requires `"T"` (temperature in K) and `"duration"` (seconds).
- `"ramp"`: requires `"rate"` (K/s); `"duration"` is optional — if omitted the ramp extends to the end of the simulation.

Example (TOML):
```toml
[temperature.tds]
type = "piecewise"
stages = [
  {type = "hold", T = 300.0, duration = 600.0},
  {type = "ramp", rate = 0.5},
]
```
"""
function PiecewiseTemperature(stages::AbstractVector)
    breakpoints = Float64[]
    T_at_bps = Float64[]
    rates = Float64[]

    t_current = 0.0
    T_current = 0.0

    for (i, stage) in enumerate(stages)
        kind = String(stage["type"])

        if kind == "hold"
            T_val = Float64(stage["T"])
            duration = Float64(stage["duration"])
            push!(breakpoints, t_current)
            push!(T_at_bps, T_val)
            push!(rates, 0.0)
            T_current = T_val
            t_current += duration
        elseif kind == "ramp"
            rate = Float64(stage["rate"])
            push!(breakpoints, t_current)
            push!(T_at_bps, T_current)
            push!(rates, rate)
            if haskey(stage, "duration")
                dur = Float64(stage["duration"])
                T_current += rate * dur
                t_current += dur
            else
                break
            end
        else
            throw(ArgumentError("Unknown temperature stage type: \"$kind\". Use \"hold\" or \"ramp\"."))
        end
    end

    isempty(breakpoints) &&
        throw(ArgumentError("PiecewiseTemperature requires at least one stage"))
    return PiecewiseTemperature(breakpoints, T_at_bps, rates)
end
