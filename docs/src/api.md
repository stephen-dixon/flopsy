# API Reference

Auto-generated from docstrings.  All symbols listed here are exported from `Flopsy`.

## Variable layout

```@docs
VariableInfo
VariableLayout
nvariables
variable_names
variables_in_group
variables_with_tag
```

## Mesh

```@docs
Mesh1D
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

## Operators

```@docs
NullOperator
OperatorSum
LinearDiffusionOperator
DirichletBoundaryOperator
ToyReactionOperator
SimpleTrappingReactionOperator
ConstraintOperator
supports_rhs
supports_jacobian
rhs!
```

## Formulations

```@docs
UnsplitFormulation
```

## Model builders

```@docs
build_rd_model
build_trapping_model
```

## Solver

```@docs
SolverConfig
build_problem
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
FakeHotgatesModel
hotgates_rates!
```

## Output

```@docs
SimulationResult
wrap_result
variable_timeseries
variable_snapshot
integrated_variable
surface_diffusive_fluxes
build_summary_dataframe
write_field_output_hdf5
write_summary_csv
load_ic_from_hdf5
```
