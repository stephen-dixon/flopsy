"""
    build_problem(model, u0, tspan, ::IMEXFormulation, solver_config) -> SplitODEProblem

Assemble a `SplitODEProblem` for IMEX integration.

- **f1 (stiff / implicit)** — diffusion and boundary operators.  When all stiff
  operators support analytic Jacobians, f1 is wrapped with a sparse `jac!` and
  prototype so that IMEX solvers can exploit the known sparsity.
- **f2 (non-stiff / explicit)** — reaction operators (and any remaining operators
  not assigned to the stiff partition).

Use with an IMEX algorithm such as `KenCarp4()`.

# Example
```julia
config = SolverConfig(
    formulation = IMEXFormulation(),
    algorithm   = KenCarp4(),
    abstol      = 1e-8,
    reltol      = 1e-6,
    saveat      = range(0.0, t_end; length = 200),
)
result = solve_problem(model, u0, (0.0, t_end), config)
```
"""
function build_problem(model::SystemModel, u0, tspan, ::IMEXFormulation, solver_config::SolverConfig)
    # Partition operators into stiff (diffusion + boundary) and non-stiff (reaction).
    stiff_ops    = Tuple(filter(!isnothing, [model.operators.diffusion,
                                             model.operators.boundary]))
    nonstiff_ops = Tuple(filter(!isnothing, [model.operators.reaction]))

    stiff_op    = isempty(stiff_ops)    ? NullOperator() : OperatorSum(stiff_ops)
    nonstiff_op = isempty(nonstiff_ops) ? NullOperator() : OperatorSum(nonstiff_ops)

    ctx = model.context

    # --- stiff part (f1) ---
    if supports_jacobian(stiff_op) && !isempty(stiff_ops)
        stiff_prototype = _build_jac_prototype(model, stiff_ops)

        function jac1!(J, u, p, t)
            fill!(J, zero(eltype(J)))
            jacobian!(J, stiff_op, u, ctx, t)
            return nothing
        end

        function f1!(du, u, p, t)
            fill!(du, zero(eltype(du)))
            rhs!(du, stiff_op, u, ctx, t)
            return nothing
        end

        f1_sciml = ODEFunction(f1!, jac = jac1!, jac_prototype = stiff_prototype)
    else
        function f1_plain!(du, u, p, t)
            fill!(du, zero(eltype(du)))
            rhs!(du, stiff_op, u, ctx, t)
            return nothing
        end
        f1_sciml = ODEFunction(f1_plain!)
    end

    # --- non-stiff part (f2) ---
    function f2!(du, u, p, t)
        fill!(du, zero(eltype(du)))
        rhs!(du, nonstiff_op, u, ctx, t)
        return nothing
    end
    f2_sciml = ODEFunction(f2!)

    return SplitODEProblem(f1_sciml, f2_sciml, u0, tspan)
end


"""
    build_problem(model, u0, tspan, ::IMEXReactionFormulation, solver_config) -> SplitODEProblem

Assemble a `SplitODEProblem` for IMEX integration with the reaction partition implicit.

- **f1 (stiff / implicit)** — reaction operators.  When the reaction operator supports
  analytic Jacobians (e.g. `HotgatesReactionOperator` backed by Palioxis), f1 is
  wrapped with a sparse `jac!` and prototype so IMEX solvers exploit the known sparsity.
- **f2 (non-stiff / explicit)** — diffusion and boundary operators.

Use with `KenCarp4()` or similar IMEX algorithms.
"""
function build_problem(model::SystemModel, u0, tspan, ::IMEXReactionFormulation,
                        solver_config::SolverConfig)
    stiff_ops    = Tuple(filter(!isnothing, [model.operators.reaction]))
    nonstiff_ops = Tuple(filter(!isnothing, [model.operators.diffusion,
                                             model.operators.boundary]))

    stiff_op    = isempty(stiff_ops)    ? NullOperator() : OperatorSum(stiff_ops)
    nonstiff_op = isempty(nonstiff_ops) ? NullOperator() : OperatorSum(nonstiff_ops)

    ctx = model.context

    # --- stiff part (f1): reaction ---
    if supports_jacobian(stiff_op) && !isempty(stiff_ops)
        stiff_prototype = _build_jac_prototype(model, stiff_ops)

        function jac1!(J, u, p, t)
            fill!(J, zero(eltype(J)))
            jacobian!(J, stiff_op, u, ctx, t)
            return nothing
        end

        function f1!(du, u, p, t)
            fill!(du, zero(eltype(du)))
            rhs!(du, stiff_op, u, ctx, t)
            return nothing
        end

        f1_sciml = ODEFunction(f1!, jac = jac1!, jac_prototype = stiff_prototype)
    else
        function f1_plain!(du, u, p, t)
            fill!(du, zero(eltype(du)))
            rhs!(du, stiff_op, u, ctx, t)
            return nothing
        end
        f1_sciml = ODEFunction(f1_plain!)
    end

    # --- non-stiff part (f2): diffusion + boundary ---
    function f2!(du, u, p, t)
        fill!(du, zero(eltype(du)))
        rhs!(du, nonstiff_op, u, ctx, t)
        return nothing
    end
    f2_sciml = ODEFunction(f2!)

    return SplitODEProblem(f1_sciml, f2_sciml, u0, tspan)
end
