# Formulations

A **formulation** determines how operators are assembled into a time-integration problem.
Flopsy provides five formulations that cover the full range from simple monolithic ODEs to
operator-split methods and DAE quasi-static approximations.

## Quick Reference

| Formulation | Type | Algorithm | When to use |
|---|---|---|---|
| `UnsplitFormulation` | Monolithic ODE | `Rodas5P()`, `CVODE_BDF()` | Default; stiff problems with analytic Jacobian |
| `IMEXFormulation` | IMEX split ODE | `KenCarp4()` | Diffusion-dominated stiffness |
| `IMEXReactionFormulation` | IMEX split ODE (reaction implicit) | `KenCarp4()` | Reaction-dominated stiffness |
| `SplitFormulation(LieSplit())` | Operator split | `Rodas5(autodiff=AutoFiniteDiff())` | Physical process separation, 1st order |
| `SplitFormulation(StrangSplit())` | Operator split | `Rodas5(autodiff=AutoFiniteDiff())` | Physical process separation, 2nd order |
| `ResidualFormulation` | Mass-matrix DAE | `Rodas5P()` | Quasi-static trapping (fast trapping limit) |

---

## UnsplitFormulation

All operators are summed into a single `f!(du, u, p, t)` and passed to a stiff ODE solver.
When all operators implement analytic Jacobians (which is the case for the built-in
`SimpleTrappingReactionOperator`, `LinearDiffusionOperator`, and `DirichletBoundaryOperator`),
`build_problem` automatically constructs a sparse `ODEFunction` with the block-tridiagonal
Jacobian structure.

For Palioxis-backed models using `HotgatesReactionOperator`, the analytic Jacobian is
unavailable and the solver falls back to finite-difference differentiation.

```julia
using Flopsy
using OrdinaryDiffEq
using ADTypes: AutoFiniteDiff

config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    abstol      = 1e-10,
    reltol      = 1e-8,
    saveat      = collect(0.0:100.0:3600.0),
)

sol    = solve_problem(model, u0, tspan, config)
result = wrap_result(model, sol, config)
```

**With Palioxis:**

```julia
using Flopsy
using Palioxis
using OrdinaryDiffEq
using ADTypes: AutoFiniteDiff

Palioxis.init(ENV["PALIOXIS_ROOT"])
pal  = MultipleDefectModel(ENV["PALIOXIS_XML"])
mesh = Mesh1D(1e-3, 200)
defects = fill(1e-3, get_n_trap_types(pal), 200)

model = build_palioxis_trapping_model(
    palioxis_model = pal,
    mesh           = mesh,
    defects        = defects,
    temperature    = LinearRampTemperature(300.0, 10.0/60.0),
    left_bc        = t -> 0.0,
    right_bc       = t -> 0.0,
)

u0 = zeros(nvariables(model.layout) * 200)

config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    abstol      = 1e-10,
    reltol      = 1e-8,
    saveat      = collect(0.0:10.0:54000.0),   # 10 s intervals over 15 h TDS ramp
)

result = wrap_result(model, solve_problem(model, u0, (0.0, 54000.0), config), config)
write_field_output_hdf5(result, "tds_fields.h5")
write_summary_csv(result, "tds_summary.csv")
```

---

## IMEXFormulation

Splits the problem into a *stiff implicit* part (diffusion + boundary operators) and a
*non-stiff explicit* part (reaction operators).  Assembles a `SplitODEProblem` and uses
an IMEX Runge-Kutta algorithm.

This can be faster than `UnsplitFormulation` when diffusion dominates stiffness and the
reaction term is relatively cheap.  It is less effective when trapping reaction timescales
are stiff too (e.g. at low temperatures with very fast trapping).

```julia
using Flopsy
using OrdinaryDiffEq

config = SolverConfig(
    formulation = IMEXFormulation(),
    algorithm   = KenCarp4(),
    abstol      = 1e-9,
    reltol      = 1e-7,
    saveat      = collect(saveat_times),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

**With Palioxis:**

The `HotgatesReactionOperator` routes to the explicit part; `LinearDiffusionOperator` and
`DirichletBoundaryOperator` form the implicit part — the partition is automatic.

```julia
config = SolverConfig(
    formulation = IMEXFormulation(),
    algorithm   = KenCarp4(),
    abstol      = 1e-10,
    reltol      = 1e-8,
    saveat      = collect(0.0:10.0:54000.0),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

---

## IMEXReactionFormulation

Like `IMEXFormulation` but with the partition swapped: the *reaction* operators form the
stiff **implicit** part and the *diffusion + boundary* operators are the **explicit** part.
Use this when trapping kinetics (not diffusion) dominate the stiffness.

When the reaction operator provides an analytic Jacobian (e.g. `SimpleTrappingReactionOperator`),
`build_problem` constructs an `ODEFunction` with a sparse Jacobian for the implicit part,
reducing the cost of each Newton iteration.

```julia
using Flopsy
using OrdinaryDiffEq

config = SolverConfig(
    formulation = IMEXReactionFormulation(),
    algorithm   = KenCarp4(),
    abstol      = 1e-9,
    reltol      = 1e-7,
    saveat      = collect(saveat_times),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

!!! note
    The implicit reaction Jacobian path is exercised when the model's `reaction` operator
    returns `true` from `supports_jacobian`.  For `HotgatesReactionOperator` this is
    `false`; the solver falls back to finite-difference differentiation for that part.

---

## SplitFormulation

Explicit operator-splitting loop with a fixed macro-step `dt`.  At each macro-step the
reaction and diffusion sub-problems are each solved as independent ODEs, then the results
are chained.  Two schemes are provided:

### LieSplit (first-order)

Each macro-step: reaction sub-step → diffusion sub-step.

```julia
using Flopsy
using OrdinaryDiffEq
using ADTypes: AutoFiniteDiff

config = SolverConfig(
    formulation = SplitFormulation(LieSplit()),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    dt          = 5.0,           # macro-step size (seconds)
    abstol      = 1e-9,
    reltol      = 1e-7,
    saveat      = collect(saveat_times),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

### StrangSplit (second-order)

Each macro-step: half reaction → full diffusion → half reaction.  Second-order accurate in
the macro-step size; preferred over `LieSplit` when a coarser `dt` is needed.

```julia
config = SolverConfig(
    formulation = SplitFormulation(StrangSplit()),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    dt          = 10.0,
    abstol      = 1e-9,
    reltol      = 1e-7,
    saveat      = collect(saveat_times),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

!!! note
    Sub-step ODEs have no analytic Jacobian.  Use
    `Rodas5(autodiff = AutoFiniteDiff())` rather than the default `Rodas5()` to avoid
    type conflicts between ForwardDiff `Dual` numbers and pre-allocated `Float64` scratch
    buffers.

**With Palioxis:**

Splitting is particularly useful when the Palioxis reaction kernel is expensive and the
diffusion step benefits from a different tolerance or step size.

```julia
using Flopsy
using Palioxis
using OrdinaryDiffEq
using ADTypes: AutoFiniteDiff

# model built as above with build_palioxis_trapping_model(...)

config = SolverConfig(
    formulation = SplitFormulation(StrangSplit()),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    dt          = 10.0,          # 10 s macro-step
    abstol      = 1e-9,
    reltol      = 1e-7,
    saveat      = collect(0.0:60.0:54000.0),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

---

## ResidualFormulation

DAE formulation using a singular mass matrix.  Variables in the `:trap` group receive
zero mass-matrix diagonals (algebraic constraint rows); all other variables are ODE rows.
This enforces *quasi-static trap equilibrium* at every time step.

Use this formulation when trapping kinetics are much faster than diffusion, so that traps
can be considered always in equilibrium with the local mobile concentration.

```julia
using Flopsy
using OrdinaryDiffEq

config = SolverConfig(
    formulation = ResidualFormulation(),
    algorithm   = Rodas5P(),
    abstol      = 1e-9,
    reltol      = 1e-7,
    saveat      = collect(saveat_times),
)

result = wrap_result(model, solve_problem(model, u0, tspan, config), config)
```

!!! warning "Different physical model"
    `ResidualFormulation` enforces quasi-static trap equilibrium.  Discrepancies between
    this and `UnsplitFormulation` are *physical* (the quasi-static approximation), not
    numerical error.  Only use it when the trap relaxation time is much shorter than the
    diffusion timescale.

---

## Choosing a Formulation

**Start with `UnsplitFormulation`.**  It uses an analytic sparse Jacobian (when available)
and is the most thoroughly tested path.

Switch formulations when:

- **`IMEXFormulation`** — you suspect diffusion-dominated stiffness and want to evaluate
  whether IMEX stepping is faster for your specific problem.
- **`SplitFormulation(StrangSplit())`** — you want to decouple physics for interpretability,
  or the Palioxis reaction kernel is expensive and you want to control the splitting step
  separately from the within-step tolerance.
- **`ResidualFormulation`** — trapping is very fast and you want to enforce the quasi-static
  equilibrium constraint to reduce problem dimensionality.

All formulations share the same model, initial condition, output helpers (`variable_timeseries`,
`build_summary_dataframe`, etc.), and HDF5/CSV writers.  You can run the same model with
multiple formulations and compare results directly.

### Example: comparing formulations on the same problem

```julia
using Flopsy, OrdinaryDiffEq, ADTypes

model   = build_trapping_model(mesh=Mesh1D(1e-3, 40), k_trap=5e-3, k_detrap=0.05,
                               diffusion_coefficient=1e-7)
u0      = # ... initial condition
tspan   = (0.0, 200.0)
saveat  = collect(range(0.0, 200.0; length=101))

configs = Dict(
    :unsplit => SolverConfig(formulation=UnsplitFormulation(),
                             algorithm=Rodas5P(), abstol=1e-9, reltol=1e-7,
                             saveat=saveat),
    :imex    => SolverConfig(formulation=IMEXFormulation(),
                             algorithm=KenCarp4(), abstol=1e-9, reltol=1e-7,
                             saveat=saveat),
    :strang  => SolverConfig(formulation=SplitFormulation(StrangSplit()),
                             algorithm=Rodas5(autodiff=AutoFiniteDiff()),
                             dt=2.0, abstol=1e-9, reltol=1e-7,
                             saveat=saveat),
)

results = Dict(name => wrap_result(model, solve_problem(model, u0, tspan, cfg), cfg)
               for (name, cfg) in configs)

# Compare final mobile profiles
for (name, r) in results
    c_end = variable_snapshot(r, :c)
    println("$name  final mobile total = $(sum(c_end) * model.context.mesh.dx)")
end
```

See `examples/formulation_comparison.jl` for a full comparison including wall-time
benchmarks and error norms.
