module Flopsy

using LinearAlgebra
using SparseArrays
using ADTypes
using SciMLBase
using OrdinaryDiffEq
using Sundials
using HDF5
using CSV
using DataFrames
using Printf
using Dates

include("core/types.jl")
include("core/variables.jl")
include("core/state.jl")
include("core/mesh.jl")
include("core/temperature.jl")
include("core/diffusion_coefficients.jl")

include("core/operators.jl")
include("core/reaction_operators.jl")
include("core/diffusion_operators.jl")
include("core/constraint_operators.jl")
include("core/operator_composition.jl")

include("core/formulations.jl")
include("core/unsplit_formulation.jl")
include("core/imex_formulation.jl")
include("core/split_formulation.jl")
include("core/residual_formulation.jl")

include("core/solver_api.jl")
include("core/problem_builders.jl")
include("core/time_integrators.jl")

include("core/models.jl")
include("core/generic_rd_model.jl")
include("core/trapping_model.jl")
include("core/hotgates_adapter.jl")

include("config/types.jl")
include("config/parse.jl")
include("config/validate.jl")
include("config/runtime_types.jl")
include("config/plugins.jl")
include("config/registry.jl")
include("config/syntax_builtins.jl")
include("config/runtime_build.jl")
include("problem_templates/diffusion_1d.jl")
include("problem_templates/trapping_1d.jl")
include("problem_templates/hotgates_trapping.jl")
include("problem_factory/build_problem.jl")

include("progress.jl")
include("output.jl")
include("plotting.jl")
include("cli.jl")

export
       AbstractSystemModel,
       AbstractOperator,
       AbstractReactionOperator,
       AbstractDiffusionOperator,
       AbstractConstraintOperator,
       AbstractFormulation, VariableInfo,
       VariableLayout,
       Mesh1D,
       SystemContext,
       SystemModel,
       SimulationProblem,
       SolverConfig,
       ParameterSpec,
       SyntaxSpec,
       SyntaxRegistry,
       BuildContext,
       ConfigBlock,
       InputDeck,
       SpeciesInfo,
       BackendDefinition,
       ICDefinition,
       BCDefinition,
       OutputDefinition,
       ProblemDefinition,
       ConfiguredSimulation,
       nvariables,
       variable_names,
       variables_in_group,
       variables_with_tag,
       has_group, state_view,
       node_view,
       variable_view,
       group_view, supports_rhs,
       supports_implicit_rhs,
       supports_step,
       supports_residual,
       supports_mass_matrix,
       supports_jacobian, rhs!,
       implicit_rhs!,
       step!,
       residual!,
       jacobian!,
       mass_matrix, NullOperator,
       OperatorSum,
       active_operators, AbstractDiffusionCoefficients,
       get_D,
       ConstantDiffusion,
       ArrheniusDiffusion,
       CallableDiffusion, LinearDiffusionOperator,
       WeakDirichletBoundaryOperator,
       DirichletBoundaryOperator,
       PenaltyMethod,
       MassMatrixMethod,
       CallbackMethod,
       EliminatedMethod,
       build_solver_callback,
       ToyReactionOperator,
       SimpleTrappingReactionOperator,
       ConstraintOperator, UnsplitFormulation,
       IMEXFormulation,
       IMEXReactionFormulation,
       LieSplit,
       StrangSplit,
       SplitFormulation,
       ResidualFormulation, build_rd_model,
       build_trapping_model,
       trapping_variable_layout, SplitProblem,
       SplitSolution, build_problem,
       solve_problem,
       solve,
       parse_input_deck,
       build_context,
       build_simulation,
       build_registry,
       register_syntax!,
       lookup_syntax,
       syntax_list,
       syntax_show,
       validate_input_deck,
       run_input_deck,
       register_plugin_provider!,
       plugin_list,
       plugin_register!,
       plugin_remove!,
       cli_main,
       julia_main,
       run_simulation, SimulationResult,
       wrap_result,
       variable_timeseries,
       variable_snapshot,
       integrated_variable,
       surface_diffusive_fluxes,
       build_summary_dataframe,
       write_field_output_hdf5,
       write_xdmf_for_flopsy_h5,
       write_summary_csv,
       load_ic_from_hdf5,
       print_run_banner,
       solver_stats_dict,
       check_mass_conservation, AbstractTemperatureProvider,
       ConstantTemperature,
       LinearRampTemperature,
       FunctionTemperature,
       temperature_at, HotgatesTrappingAdaptor,
       HotgatesReactionOperator,
       build_hotgates_variable_layout,
       build_hotgates_trapping_model,
       build_palioxis_trapping_model,
       build_equilibrium_ic,
       build_ic_from_total_hydrogen, FakeHotgatesModel,
       hotgates_rates!, plot_tds_flux,
       plot_spatial_snapshot,
       plot_spatial_evolution,
       record_spatial_video,
       record_spatial_animation
end
