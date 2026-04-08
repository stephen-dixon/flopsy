# Architecture

## State Vector Layout

Flopsy models evolve a flat state vector `u` of length `nvars * nx`, where `nvars` is
the number of physical variables and `nx` is the number of spatial nodes.

The vector is packed in **variable-major (column-major) order**: all variables at node 1
come first, then all variables at node 2, and so on:

```
u[1]              ← Variable 1, node 1
u[2]              ← Variable 2, node 1
⋮
u[nvars]          ← Variable nvars, node 1
u[nvars + 1]      ← Variable 1, node 2
⋮
u[2*nvars]        ← Variable nvars, node 2
⋮
```

The flat index of variable `ivar` at node `ix` is `(ix-1)*nvars + ivar`.

A `VariableLayout` describes the structure: each variable has a **name** (e.g. `:mobile_H`),
a **group** (e.g. `:mobile`, `:trap`), and **tags** (e.g. `:reaction`, `:diffusion`) that
control which operators act on it.

### Reshaping and Views

`state_view(u, layout, nx)` reshapes the flat vector into a `(nvars, nx)` matrix without
copying, so `U[ivar, ix]` accesses variable `ivar` at node `ix`.

Further helpers work on the already-reshaped matrix `U`:

- `node_view(U, ix)` — column `ix`; the full state at one node.
- `variable_view(U, ivar)` — row `ivar`; one variable across all nodes.
- `group_view(U, layout, group)` — contiguous rows belonging to a named group.

These views are zero-copy and are used pervasively inside operators and the Palioxis bridge.

## Operators

Flopsy defines a hierarchy of abstract operator types:

- `AbstractReactionOperator` — node-local reaction rates (trapping, detrapping, decay).
- `AbstractDiffusionOperator` — spatial transport and boundary conditions.
- `AbstractConstraintOperator` — algebraic constraints (DAE path, currently scaffolded).

### Capability Flags

Each operator declares what it implements via predicate functions:

| Predicate | When `true` |
|---|---|
| `supports_rhs(op)` | `rhs!(du, op, u, ctx, t)` accumulates an explicit RHS contribution |
| `supports_jacobian(op)` | `jacobian!(J, op, u, ctx, t)` accumulates analytic ∂f/∂u entries |
| `supports_implicit_rhs(op)` | operator contributes to an implicit split |
| `supports_mass_matrix(op)` | operator provides a mass matrix (DAE path) |

Defaults are `false`; operators opt in by adding the relevant method.

### Analytic Jacobians

`LinearDiffusionOperator`, `DirichletBoundaryOperator`, `ToyReactionOperator`, and
`SimpleTrappingReactionOperator` all implement `supports_jacobian = true` with analytic
`jacobian!` methods.

The Jacobian structure exploits the known sparsity:

- **Diffusion**: tridiagonal coupling along each diffusing variable's spatial dimension
  (same variable, adjacent nodes).
- **Reaction**: dense within each node's variable block (all variables at a node can
  couple to each other), but no inter-node coupling.
- **Combined**: block-tridiagonal sparse matrix, built automatically as a `SparseArrays`
  prototype by `build_unsplit_problem` when `supports_jacobian` returns `true` for the
  assembled `OperatorSum`.

For models where Jacobian support is unavailable (e.g. `HotgatesReactionOperator`),
`build_unsplit_problem` falls back to `ODEFunction(f!)` without an explicit Jacobian,
relying on the algorithm's own differentiation strategy.

### The RHS Contract

```julia
rhs!(du, op, u, ctx, t)
```

- `du` — pre-allocated output vector; the operator **accumulates** into it (does not zero it).
- `u`  — current flat state vector.
- `ctx` — `SystemContext` with mesh, layout, and scratch buffers.
- `t`  — current time.

`OperatorSum.rhs!` zeros `du` before calling each sub-operator, so the total RHS is the
sum of all operator contributions.

### DirichletBoundaryOperator

Boundary conditions are operators. `DirichletBoundaryOperator` adds ghost-node Dirichlet
corrections on top of the zero-flux stencil already in `LinearDiffusionOperator`:

```julia
boundary = DirichletBoundaryOperator(
    selector, coefficients, temperature;
    left  = t -> 0.0,
    right = t -> 0.0,
)
```

The combined stencil at node 1 becomes `D*(U[2] - 2*U[1] + g)/dx²`, where `g = left(t)`.

### Composition with OperatorSum

Complex problems combine operators via `OperatorSum`:

```julia
total = OperatorSum((reaction_op, diffusion_op, boundary_op))
rhs!(du, total, u, ctx, t)    # calls each sub-operator and accumulates
jacobian!(J, total, u, ctx, t) # sums analytic Jacobian contributions (if all support it)
```

`active_operators(model)` returns all non-`nothing` slots from `model.operators` as a
`Vector`, which `build_unsplit_problem` wraps in a `OperatorSum`.

## Diffusion Coefficients

`LinearDiffusionOperator` and `DirichletBoundaryOperator` accept either:

- A plain `Vector{<:Real}` — constant, temperature-independent.
- Any `AbstractDiffusionCoefficients` subtype — evaluated via `get_D(coeffs, ivar, ix, T)`:
  - `ConstantDiffusion` — ignores both `ix` and `T`.
  - `ArrheniusDiffusion` — D(T) = D₀ exp(−Eₐ/kBT), ignores `ix`.
  - `CallableDiffusion` — one callable `f(T)` per variable.
  - `PalioxisDiffusionCoefficients` (Palioxis extension) — queries `Palioxis.diffusion_constants` at every evaluation.

The `ix` argument is reserved for future spatially-varying D; all current types ignore it.

## Formulations

A **formulation** specifies how operators are assembled into a SciML-compatible problem.

### UnsplitFormulation

The default: all active operators are summed into a single ODE.

```
du/dt = Σ rhs!(op, u, t)   for all active operators
```

When `supports_jacobian` is true for the assembled `OperatorSum`, `build_unsplit_problem`
automatically constructs a sparse `ODEFunction` with an analytic Jacobian:

```julia
ODEFunction(f!, jac = jac!, jac_prototype = sparse_prototype)
```

The sparse prototype is built from the known structure: dense within-node blocks (reaction
coupling) plus tridiagonal off-diagonals (diffusion coupling).

### Other Formulations (Stubs)

- `IMEXFormulation` — separate explicit and implicit parts for semi-implicit integration.
- `SplitFormulation` / `LieSplit` / `StrangSplit` — operator-split time stepping.
- `ResidualFormulation` — residual-based DAE path.

These are scaffolded and will be implemented in a future pass.  Production simulations
should use `UnsplitFormulation`.

## SystemModel and SystemContext

### SystemModel

Bundles the complete problem description:

- **VariableLayout** — structure of the state vector.
- **Operators** — named tuple `(reaction, diffusion, boundary, constraints)`.  Each slot
  may be `nothing`.
- **SystemContext** — mesh, node count, aux data, scratch buffers.

Built via `build_rd_model`:

```julia
model = build_rd_model(
    layout      = layout,
    mesh        = mesh,
    reaction    = reaction_op,
    diffusion   = diffusion_op,
    boundary    = boundary_op,   # DirichletBoundaryOperator or nothing
)
```

### SystemContext

Holds:

- **`layout`** — `VariableLayout`.
- **`nx`** — number of spatial nodes.
- **`mesh`** — `Mesh1D` with coordinates and spacing.
- **`aux`** — `Dict{Symbol,Any}` for auxiliary data.
- **`scratch`** — `Dict{Symbol,Any}` for pre-allocated working arrays (populated lazily
  by operators via `get!(ctx.scratch, key) do ... end`).

## SolverConfig and Solve Pipeline

`SolverConfig` wraps:

- **`formulation`** — typically `UnsplitFormulation()`.
- **`algorithm`** — SciML algorithm, e.g. `Rodas5(autodiff=AutoFiniteDiff())`.
- **`abstol`, `reltol`** — solver tolerances.
- **`saveat`** — output times.

### Pipeline

```julia
sol    = solve_problem(model, u0, tspan, solver_config)
result = wrap_result(model, sol, config)
```

`solve_problem` calls `build_problem`, which dispatches on the formulation to construct
the `ODEProblem` (with sparse Jacobian if available), then calls `SciMLBase.solve`.
