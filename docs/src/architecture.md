# Architecture

## State Vector Layout

Flopsy models evolve a flat state vector `u` of length `nvars * nx`, where `nvars` is the number of physical variables and `nx` is the number of spatial nodes. A `VariableLayout` describes how each variable is packed into this vector:

```
u[1:nx]           ← Variable 1, nodes 1:nx
u[nx+1:2*nx]      ← Variable 2, nodes 1:nx
...
```

Each variable in the layout has:
- A **name** (e.g., `"mobile_H"`, `"trapped_H_defect1"`).
- A **group** (`:mobile` or `:trap`) — indicates whether the variable represents free or trapped hydrogen.
- **Tags** (e.g., `:reaction`, `:diffusion`) — indicate which operators act on this variable.

### Reshaping and Views

The helper function `state_view(u, layout, nx)` reshapes `u` into a `(nvars, nx)` matrix, making it easier to index by variable and node. Further helper functions provide intuitive access:

- `node_view(u, layout, ix)` — extract the state at a single spatial node `ix`.
- `variable_view(u, layout, name)` — extract all values of a variable across all nodes.
- `group_view(u, layout, group)` — extract all variables in a given group (e.g., all trapped species).

These views are essential for localizing computations and for passing state subsets to external libraries like Palioxis.

## Operators

The package provides a set of abstract operator types that encapsulate different physics:

- `AbstractReactionOperator` — governs reaction rates (trapping, detrapping, etc.).
- `AbstractDiffusionOperator` — governs spatial diffusion and boundary conditions.
- `AbstractConstraintOperator` — enforces algebraic constraints (e.g., conservation laws, flux boundary conditions).

### Capability Flags

Each operator declares what it can compute via capability flags:

- `supports_rhs(op)` — whether `rhs!(du, op, u, ctx, t)` is implemented.
- `supports_jacobian(op)` — whether Jacobian blocks can be computed (for implicit solvers).
- `supports_mass_matrix(op)` — whether a mass-matrix is defined (for DAE systems).

### The RHS Contract

All operators implement a common interface:

```julia
rhs!(du, op, u, ctx, t)
```

where:
- `du` is the rate-of-change vector (in-place update).
- `op` is the operator instance.
- `u` is the current state vector.
- `ctx` is a `SystemContext` holding mesh, auxiliary data, and scratch space.
- `t` is the current time.

The operator **accumulates** contributions to `du` (does not overwrite).

### Example: DirichletBoundaryOperator

Boundary conditions are themselves operators. The `DirichletBoundaryOperator` adds time-varying Dirichlet corrections on top of the zero-flux stencil:

```julia
boundary = DirichletBoundaryOperator(
    selector,
    coefficients;
    left = t -> 0.0,      # callable boundary value
    right = t -> 0.0,
)
```

The operator queries the callables `left(t)` and `right(t)` at each time step and adds the corrections:

```
dU[ivar, 1]  += D * (left(t)  - U[ivar, 1])  / dx²
dU[ivar, nx] += D * (right(t) - U[ivar, nx]) / dx²
```

Combined with the centred stencil from `LinearDiffusionOperator`, this gives the standard Dirichlet discretisation.

### Composition with OperatorSum

Complex problems combine multiple operators via `OperatorSum`:

```julia
op = OperatorSum([reaction_op, diffusion_op, boundary_op, constraint_op])
rhs!(du, op, u, ctx, t)  # calls rhs! on each operator and sums the results
```

Each operator accumulates its contribution independently, making the composition modular and extensible.

## Formulations

A **formulation** specifies how operators are combined to construct a SciML-compatible problem. Different formulations support different time-integration strategies:

### UnsplitFormulation

The production formulation: all operators are combined into a single ODE system

```
du/dt = reaction(u, t) + diffusion(u, t) + boundary(u, t) + ...
```

The operators are summed directly in `rhs!`, and the entire state is integrated implicitly. This fully-implicit approach is robust for stiff systems.

### Other Formulations (Stubs)

- `IMEXFormulation` — separate explicit and implicit operators for semi-implicit time integration.
- `OperatorSplitFormulation` — time-split the operators (e.g., reaction steps followed by diffusion steps).
- `DAEFormulation` — residual-based formulation for systems with algebraic constraints.

These are currently scaffolded; production use should employ `UnsplitFormulation`.

## SystemModel and SystemContext

### SystemModel

A `SystemModel` bundles the complete problem description:

- **VariableLayout** — the structure of the state vector.
- **Operators** — a named tuple `(reaction, diffusion, boundary, constraints)`. Each slot may be `nothing` if that physics is absent. The `boundary` slot was added to support operators that modify boundary behavior independently of the main diffusion stencil.
- **SystemContext** — mesh, spatial step size, auxiliary parameters, and scratch arrays.

It serves as the central object passed through the solve pipeline.

```julia
model = SystemModel(layout, operators, context)
```

The operators tuple is created via `build_rd_model`:

```julia
model = build_rd_model(
    layout = layout,
    mesh = mesh,
    reaction = reaction_op,
    diffusion = diffusion_op,
    boundary = boundary_op,        # DirichletBoundaryOperator or nothing
    constraints = constraint_op,
    aux = Dict(:temp_profile => ..., ...)
)
```

### SystemContext

The `SystemContext` holds:

- **layout** — `VariableLayout` describing state vector structure.
- **nx** — number of spatial nodes.
- **mesh** — `Mesh1D` with grid coordinates and spacing.
- **aux** — dictionary of auxiliary data (temperature profiles, material properties, defect densities).
- **scratch** — pre-allocated arrays for in-place computations (populated as needed by operators).

By centralizing auxiliary data and pre-allocated arrays, `SystemContext` ensures that operators and time integrators can compute efficiently without repeated allocations.

## SolverConfig and Solve Pipeline

A `SolverConfig` wraps:

- **Formulation** — which formulation to use (typically `UnsplitFormulation`).
- **Algorithm** — the time-integration algorithm (e.g., `Rodas5`, `CVODE_BDF` via Sundials.jl).
- **Tolerances** — absolute and relative error tolerances.
- **Save settings** — output frequency and dense output options.

### The Solve Pipeline

1. **build_problem** — Given a `SystemModel` and `SolverConfig`, construct a SciML `ODEProblem`.
2. **solve_problem** — Solve the problem using the specified algorithm and tolerances.

The pipeline abstracts away SciML boilerplate: the user provides a high-level model description, and Flopsy handles the conversion to a form suitable for robust implicit solvers.

```julia
model = SystemModel(layout, operators, context)
config = SolverConfig(formulation, algorithm, rtol, atol)
problem = build_problem(model, config)
solution = solve_problem(problem, config)
```

This design ensures that adding new operators or formulations requires minimal changes to the solver infrastructure.
