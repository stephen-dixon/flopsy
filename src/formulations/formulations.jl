"""
    UnsplitFormulation

Monolithic ODE formulation: all operators are summed into a single `f!(du,u,p,t)`
and passed to a stiff ODE solver (e.g. Rodas5, CVODE_BDF).  This is the primary
production path in Flopsy.
"""
struct UnsplitFormulation <: AbstractFormulation end

"""
    IMEXFormulation

IMEX (implicit-explicit) splitting formulation.  Stiff terms are treated implicitly;
non-stiff terms explicitly.  *Not yet implemented — stub only.*
"""
struct IMEXFormulation <: AbstractFormulation end

"""
    LieSplit

First-order (Lie) operator-splitting scheme for use with `SplitFormulation`.
*Not yet implemented — stub only.*
"""
struct LieSplit end

"""
    StrangSplit

Second-order (Strang) operator-splitting scheme for use with `SplitFormulation`.
*Not yet implemented — stub only.*
"""
struct StrangSplit end

"""
    SplitFormulation(scheme)

Operator-splitting formulation parameterised by a splitting scheme (`LieSplit` or
`StrangSplit`).  *Not yet implemented — stub only.*
"""
struct SplitFormulation{S} <: AbstractFormulation
    scheme::S
end

"""
    ResidualFormulation

DAE residual formulation for algebraic constraints (e.g. Dirichlet BCs via
`ConstraintOperator`).  *Not yet implemented — stub only.*
"""
struct ResidualFormulation <: AbstractFormulation end
