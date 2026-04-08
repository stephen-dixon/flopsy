"""
    build_problem(model, u0, tspan, ::ResidualFormulation, solver_config) -> ODEProblem

Assemble an `ODEProblem` with a **singular mass matrix** that treats trap variables
as algebraic (quasi-static) rather than dynamic.

Variables in the `:trap` group are assigned mass-matrix diagonal entry **0**,
meaning their `du/dt` row is interpreted as an algebraic constraint `0 = f(u)`.
All other variables have diagonal entry **1** (ordinary ODEs).

This formulation is appropriate when trapping is fast relative to diffusion, so
trapped species are always at quasi-static equilibrium with the local mobile
concentration.  Use a DAE-capable algorithm such as `Rodas5P()`.

!!! note
    The RHS function `f!` is the same as in `UnsplitFormulation`.  The solver
    interprets the mass-matrix zero rows as constraints `0 = f_trap(u)`, so the
    trap reaction operator must supply `dθ/dt = 0` residuals that encode the
    equilibrium condition.  With `SimpleTrappingReactionOperator` this means the
    reaction rates will be driven to zero at steady state, not that the rates are
    literally zero — the DAE solver enforces consistency.

# Example
```julia
config = SolverConfig(
    formulation = ResidualFormulation(),
    algorithm   = Rodas5P(),
    abstol      = 1e-8,
    reltol      = 1e-6,
    saveat      = range(0.0, t_end; length = 200),
)
result = solve_problem(model, u0, (0.0, t_end), config)
```
"""
function build_problem(model::SystemModel, u0, tspan, ::ResidualFormulation, solver_config::SolverConfig)
    ops      = Tuple(active_operators(model))
    total_op = OperatorSum(ops)

    supports_rhs(total_op) || throw(ArgumentError(
        "ResidualFormulation requires rhs! support for all active operators"
    ))

    ctx = model.context

    function f!(du, u, p, t)
        rhs!(du, total_op, u, ctx, t)
        return nothing
    end

    M = _build_mass_matrix(model)

    if supports_jacobian(total_op)
        prototype = _build_jac_prototype(model, ops)

        function jac!(J, u, p, t)
            jacobian!(J, total_op, u, ctx, t)
            return nothing
        end

        ode_f = ODEFunction(f!, jac = jac!, jac_prototype = prototype, mass_matrix = M)
    else
        ode_f = ODEFunction(f!, mass_matrix = M)
    end

    return ODEProblem(ode_f, u0, tspan)
end


"""
    _build_mass_matrix(model) -> Diagonal

Build a diagonal mass matrix for the state vector:
- Diagonal entry = **1** for variables outside the `:trap` group (standard ODE).
- Diagonal entry = **0** for variables in the `:trap` group (algebraic constraint).
"""
function _build_mass_matrix(model::SystemModel)
    layout = model.context.layout
    nx     = model.context.nx
    nvars  = nvariables(layout)
    n      = nvars * nx

    # Identify which variable indices are in the :trap group.
    trap_ivars = Set(variables_in_group(layout, :trap))

    diag = ones(Float64, n)
    for ix in 1:nx
        offset = (ix - 1) * nvars
        for ivar in trap_ivars
            diag[offset + ivar] = 0.0
        end
    end

    return Diagonal(diag)
end
