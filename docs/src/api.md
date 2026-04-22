# API Reference

Auto-generated from docstrings.

## Stable Public API

### Core abstract types

```@docs
AbstractSystemModel
AbstractOperator
AbstractReactionOperator
AbstractDiffusionOperator
AbstractConstraintOperator
AbstractFormulation
```

### Variable layout

```@docs
VariableInfo
VariableLayout
nvariables
variable_names
has_group
variables_in_group
variables_with_tag
```

### Mesh and system model

```@docs
Mesh1D
SystemContext
SystemModel
SimulationProblem
SolverConfig
```

### Registry-driven input deck

```@docs
ParameterSpec
SyntaxSpec
SyntaxRegistry
BuildContext
ConfigBlock
InputDeck
BackendDefinition
ICDefinition
BCDefinition
OutputDefinition
ProblemDefinition
ConfiguredSimulation
parse_input_deck
build_registry
register_syntax!
lookup_syntax
syntax_list
syntax_show
build_context
build_simulation
validate_input_deck
run_input_deck
register_plugin_provider!
plugin_list
plugin_register!
plugin_remove!
cli_main
julia_main
```

### State vector helpers

```@docs
state_view
node_view
variable_view
group_view
```

### Temperature providers

```@docs
AbstractTemperatureProvider
ConstantTemperature
LinearRampTemperature
FunctionTemperature
PiecewiseTemperature
temperature_at
```

### Diffusion coefficients

```@docs
AbstractDiffusionCoefficients
get_D
ConstantDiffusion
ArrheniusDiffusion
CallableDiffusion
```

### Operators

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

### Formulations

```@docs
UnsplitFormulation
IMEXFormulation
IMEXReactionFormulation
LieSplit
StrangSplit
SplitFormulation
ResidualFormulation
```

### Model builders and solver core

```@docs
build_rd_model
build_trapping_model
trapping_variable_layout
SplitProblem
SplitSolution
build_problem(::SystemModel, Any, Any, ::UnsplitFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::IMEXFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::IMEXReactionFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::SplitFormulation, ::SolverConfig)
build_problem(::SystemModel, Any, Any, ::ResidualFormulation, ::SolverConfig)
solve_problem
solve(::SimulationProblem)
run_simulation
```

### Hotgates adapter

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

### Output and plotting

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
plot_tds_flux
plot_spatial_snapshot
plot_spatial_evolution
record_spatial_video
record_spatial_animation
```

## Deprecated Legacy API

These symbols still exist for compatibility but are deprecated and are no longer part of the recommended TOML workflow:

```@docs
Flopsy.MeshConfig
Flopsy.BoundaryConditionConfig
Flopsy.InputSolverConfig
Flopsy.InitialConditionConfig
Flopsy.ProblemConfig
Flopsy.load_config
Flopsy.parse_config
Flopsy.validate
Flopsy.build_problem(::Flopsy.ProblemConfig)
```
