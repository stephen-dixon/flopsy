"""
    UnsplitFormulation

Monolithic ODE formulation: all operators are summed into a single `f!(du,u,p,t)`
and passed to a stiff ODE solver (e.g. Rodas5, CVODE_BDF).  This is the primary
production path in Flopsy.
"""
struct UnsplitFormulation <: AbstractFormulation end

"""
    IMEXFormulation

IMEX (implicit-explicit) splitting formulation.  Diffusion and boundary operators
are treated implicitly (stiff partition); reaction operators are treated explicitly
(non-stiff partition).  Assembles a `SplitODEProblem` for use with IMEX algorithms
such as `KenCarp4()`.
"""
struct IMEXFormulation <: AbstractFormulation end

"""
    LieSplit

First-order (Lie) operator-splitting scheme for use with `SplitFormulation`.
Each macro-step applies: reaction sub-step → diffusion sub-step.
"""
struct LieSplit end

"""
    StrangSplit

Second-order (Strang) operator-splitting scheme for use with `SplitFormulation`.
Each macro-step applies: half reaction → full diffusion → half reaction.
"""
struct StrangSplit end

"""
    SplitFormulation(scheme)

Operator-splitting formulation parameterised by a splitting scheme (`LieSplit` or
`StrangSplit`).  Requires `solver_config.dt` to be set (macro-step size).

Sub-steps are solved independently using the configured algorithm.  Because
sub-step `ODEFunction`s have no analytic Jacobian, use an algorithm compatible
with finite-difference Jacobians, e.g. `Rodas5(autodiff=AutoFiniteDiff())`.
"""
struct SplitFormulation{S} <: AbstractFormulation
    scheme::S
end

"""
    ResidualFormulation

DAE formulation using a singular mass matrix.  Variables in the `:trap` group are
assigned a zero mass-matrix diagonal (algebraic constraint row); all other variables
are standard ODE rows (diagonal = 1).

This enforces quasi-static trap equilibrium and is appropriate when trapping kinetics
are fast relative to diffusion.  Use a DAE-capable algorithm such as `Rodas5P()`.

!!! note
    This is a *different physical model* from `UnsplitFormulation`.  Differences in
    solution reflect the quasi-static approximation, not solver error.
"""
struct ResidualFormulation <: AbstractFormulation end
