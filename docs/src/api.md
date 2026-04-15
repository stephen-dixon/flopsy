# API Reference

Auto-generated from docstrings.  All symbols listed here are exported from `Flopsy`.

## Core abstract types

```@docs
AbstractSystemModel
AbstractOperator
AbstractReactionOperator
AbstractDiffusionOperator
AbstractConstraintOperator
AbstractFormulation
```

## Variable layout

```@docs
VariableInfo
VariableLayout
nvariables
variable_names
has_group
variables_in_group
variables_with_tag
```

## Mesh

```@docs
Mesh1D
```

## System model

```@docs
SystemContext
SystemModel
SimulationProblem
```

## Config layer

```@docs
MeshConfig
BoundaryConditionConfig
InputSolverConfig
InitialConditionConfig
ProblemConfig
load_config
parse_config
validate
build_problem(::ProblemConfig)
solve(::SimulationProblem)
run_simulation
```

## State vector helpers

```@docs
state_view
node_view
variable_view
group_view
```

## Temperature providers

```@docs
AbstractTemperatureProvider
ConstantTemperature
LinearRampTemperature
FunctionTemperature
temperature_at
```

## Diffusion coefficients

```@docs
AbstractDiffusionCoefficients
get_D
ConstantDiffusion
ArrheniusDiffusion
CallableDiffusion
```

## Operators

```@docs
NullOperator
OperatorSum
active_operators
LinearDiffusionOperator
WeakDirichletBoundaryOperator
DirichletBoundaryOperator
PenaltyMethod
MassMatrixMethod
CallbackMethod
EliminatedMethod
build_solver_callback
ToyReactionOperator
SimpleTrappingReactionOperator
ConstraintOperator
supports_rhs
supports_implicit_rhs
supports_step
supports_residual
supports_mass_matrix
supports_jacobian
rhs!
implicit_rhs!
step!
residual!
jacobian!
mass_matrix
```

## Formulations

```@docs
UnsplitFormulation
IMEXFormulation
IMEXReactionFormulation
LieSplit
StrangSplit
SplitFormulation
ResidualFormulation
```

## Model builders

```@docs
build_rd_model
build_trapping_model
trapping_variable_layout
```

## Solver core

```@docs
SolverConfig
SplitProblem
SplitSolution
build_problem(::SystemModel, Any, Any, ::UnsplitFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::IMEXFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::IMEXReactionFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::SplitFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::ResidualFormulation, ::SolverConfig)
solve_problem
```

## Hotgates adapter

```@docs
HotgatesTrappingAdaptor
HotgatesReactionOperator
build_hotgates_variable_layout
build_hotgates_trapping_model
build_palioxis_trapping_model
build_equilibrium_ic
build_ic_from_total_hydrogen
FakeHotgatesModel
hotgates_rates!
```

## Plotting

```@docs
plot_tds_flux
plot_spatial_snapshot
plot_spatial_evolution
record_spatial_video
record_spatial_animation
```

## Output

```@docs
SimulationResult
wrap_result
variable_timeseries
variable_snapshot
integrated_variable
surface_diffusive_fluxes
check_mass_conservation
build_summary_dataframe
write_field_output_hdf5
write_xdmf_for_flopsy_h5
write_summary_csv
load_ic_from_hdf5
print_run_banner
solver_stats_dict
```
