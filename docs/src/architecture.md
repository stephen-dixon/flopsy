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

`LinearDiffusionOperator`, `WeakDirichletBoundaryOperator`, `ToyReactionOperator`, and
`SimpleTrappingReactionOperator` all implement `supports_jacobian = true` with analytic
`jacobian!` methods.  `DirichletBoundaryOperator` provides analytic Jacobians for all
methods except `EliminatedMethod` (whose cross-node entries are handled symbolically).

The Jacobian structure exploits the known sparsity:

- **Diffusion**: tridiagonal coupling along each diffusing variable's spatial dimension
  (same variable, adjacent nodes).
- **Reaction**: dense within each node's variable block (all variables at a node can
  couple to each other), but no inter-node coupling.
- **Combined**: block-tridiagonal sparse matrix, built automatically as a `SparseArrays`
  prototype by `build_unsplit_problem` when `supports_jacobian` returns `true` for the
  assembled `OperatorSum`.

#### Per-node sparsity patterns (`jacobian_node_sparsity`)

Each operator can declare its *per-node* Jacobian non-zero pattern via:

```julia
jacobian_node_sparsity(op, layout) -> Set{Tuple{Int,Int}} or nothing
```

The returned set contains `(row_var, col_var)` index pairs for entries that can be
non-zero within a single node's variable block.  Returning `nothing` signals that the
operator's per-node pattern is dense (or unknown) and forces a dense fallback.

`_build_jac_prototype` unions the patterns from all active operators; if any returns
`nothing`, the full dense within-node block is used.  This mechanism is used by
`HotgatesReactionOperator` to expose a selective trap–trap tridiagonal sparsity when
defect groups are provided through the `trap_groups` field of `HotgatesTrappingAdaptor`.

For models where Jacobian support is unavailable (e.g. `HotgatesReactionOperator`
without a `trap_groups` structure), `build_unsplit_problem` falls back to `ODEFunction(f!)`
without an explicit Jacobian, relying on the algorithm's own differentiation strategy.

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

### Boundary Condition Operators

#### WeakDirichletBoundaryOperator (ghost-node, recommended default)

Adds ghost-node corrections on top of the zero-flux stencil already in
`LinearDiffusionOperator`:

```julia
boundary = WeakDirichletBoundaryOperator(selector, coefficients, temperature;
                                          left  = t -> 0.0,
                                          right = t -> 0.0)
```

The combined stencil at node 1 becomes `D*(U[2] - 2*U[1] + g)/dx²`.  Both sides are
controlled from a single operator; pass `left=nothing` or `right=nothing` to leave a
surface at zero-flux (Neumann).

#### DirichletBoundaryOperator (strong, per-side)

Strong Dirichlet condition applied to one side at a time.  Four enforcement strategies
are available (set via the `method` keyword):

| Method | Mechanism | Notes |
|---|---|---|
| `PenaltyMethod(λ)` | `du[bdy] += λ*(g-u)` | Adds stiffness ∝ λ; compatible with all solvers |
| `MassMatrixMethod()` | Algebraic constraint (DAE) | Requires `Rodas5P()` or similar |
| `CallbackMethod()` | Discrete callback resets `u[bdy]` each step | Exact at step boundaries |
| `EliminatedMethod()` | Cancels Neumann stencil; sets `du[bdy]=g'(t)` | Best for constant/slow BCs |

```julia
# Left surface: strong penalty; right surface: ghost-node weak
left_bc  = DirichletBoundaryOperator(:left,  t -> 0.0, selector; method=PenaltyMethod(1e12))
right_bc = WeakDirichletBoundaryOperator(selector, coefficients; right = t -> 0.0)
boundary = OperatorSum((left_bc, right_bc))
```

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

### IMEX splitting

Partitions the operator set into a *stiff implicit* part (diffusion + boundary operators)
and a *non-stiff explicit* part (reaction operators).  Assembles a `SplitODEProblem` for
use with IMEX algorithms such as `KenCarp4()`.

Useful when the stiff timescales come primarily from diffusion and the reaction step is
cheap to handle explicitly.

```julia
config = SolverConfig(
    formulation = IMEXFormulation(),
    algorithm   = KenCarp4(),
    abstol = 1e-9, reltol = 1e-7,
)
```

### Operator splitting

Explicit operator-splitting loop with a fixed macro-step size `dt`.  Two schemes are
available:

- **`LieSplit`** (first-order): reaction → diffusion per macro-step.
- **`StrangSplit`** (second-order): half reaction → full diffusion → half reaction.

Each sub-step is solved independently as its own ODE using the configured algorithm.
Because sub-step `ODEFunction`s have no analytic Jacobian, use an algorithm that accepts
finite-difference differentiation, e.g. `Rodas5(autodiff=AutoFiniteDiff())`.

`solver_config.dt` **must** be set; it controls the splitting macro-step.

```julia
config = SolverConfig(
    formulation = SplitFormulation(StrangSplit()),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    dt          = 2.0,           # macro-step size (seconds)
    abstol = 1e-9, reltol = 1e-7,
)
```

`build_problem` returns a `SplitProblem`; `solve_problem` executes the splitting loop and
returns a `SplitSolution`.  Both expose `.t` and `.u` fields so they are compatible with
all Flopsy output helpers.

### DAE / residual formulation

DAE path using a singular mass matrix.  Variables in the `:trap` group receive a zero
diagonal in the mass matrix (algebraic rows), enforcing quasi-static trap equilibrium.
All other variables are standard ODE rows.

Use a DAE-capable algorithm such as `Rodas5P()`.

```julia
config = SolverConfig(
    formulation = ResidualFormulation(),
    algorithm   = Rodas5P(),
    abstol = 1e-9, reltol = 1e-7,
)
```

!!! note
    `ResidualFormulation` is a *different physical model* from `UnsplitFormulation`.
    Differences in the solution reflect the quasi-static approximation, not solver error.
    It is appropriate only when trapping kinetics are fast relative to diffusion.

See **[Formulations](formulations.md)** for a full comparison with worked examples.

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
- **`aux`** — auxiliary model data.
- **`scratch`** — pre-allocated working arrays and workspaces used by operators.

In the refactored codebase, `scratch` is part of the core performance boundary: operator
composition and the Hotgates bridge reuse buffers from `ctx.scratch` instead of
allocating new arrays in hot paths.

## SolverConfig and Solve Pipeline

`SolverConfig` wraps:

- **`formulation`** — typically `UnsplitFormulation()`.
- **`algorithm`** — SciML algorithm, e.g. `Rodas5(autodiff=AutoFiniteDiff())`.
- **`abstol`, `reltol`** — solver tolerances.
- **`saveat`** — output times.
- **`dt`** — optional fixed time step or split macro-step.

### Pipeline

```julia
sol    = solve_problem(model, u0, tspan, solver_config)
result = wrap_result(model, sol, config)
```

`solve_problem` calls `build_problem`, which dispatches on the formulation to construct
the `ODEProblem` (with sparse Jacobian if available), then calls `SciMLBase.solve`.

## Config Layer Boundary

The core architecture above is intentionally separate from the config-driven problem
layer:

- `load_config(path)` parses TOML into typed config structs
- `validate(cfg)` checks the config before solver construction
- `build_problem(cfg::ProblemConfig)` assembles a `SimulationProblem`
- `solve(problem::SimulationProblem)` executes the core solve path

This keeps the solver framework extensible for direct Julia use while still supporting
input-deck workflows for built-in problem templates.
