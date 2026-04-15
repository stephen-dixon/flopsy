module Flopsy

using LinearAlgebra
using SparseArrays
using ADTypes
using TOML
using SciMLBase
using OrdinaryDiffEq
using HDF5
using CSV
using DataFrames
using Printf
using Dates

include("types.jl")
include("variables.jl")
include("state.jl")
include("mesh.jl")
include("config.jl")
include("temperature.jl")
include("diffusion_coefficients.jl")

include("operators/operators.jl")
include("operators/reaction.jl")
include("operators/diffusion.jl")
include("operators/constraints.jl")
include("operators/composition.jl")

include("formulations/formulations.jl")
include("formulations/unsplit.jl")
include("formulations/imex.jl")
include("formulations/split.jl")
include("formulations/residual.jl")

include("solvers/solvers.jl")
include("solvers/sciML_problem_builders.jl")
include("solvers/time_integrators.jl")

include("models/models.jl")
include("models/generic_rd_model.jl")
include("models/trapping_model.jl")

include("adapters/hotgates.jl")

include("progress.jl")
include("runner.jl")
include("output.jl")
include("plotting.jl")

export
    AbstractSystemModel,
    AbstractOperator,
    AbstractReactionOperator,
    AbstractDiffusionOperator,
    AbstractConstraintOperator,
    AbstractFormulation,

    VariableInfo,
    VariableLayout,
    Mesh1D,
    SystemContext,
    SystemModel,
    SolverConfig,

    nvariables,
    variable_names,
    variables_in_group,
    variables_with_tag,
    has_group,

    state_view,
    node_view,
    variable_view,
    group_view,

    supports_rhs,
    supports_implicit_rhs,
    supports_step,
    supports_residual,
    supports_mass_matrix,
    supports_jacobian,

    rhs!,
    implicit_rhs!,
    step!,
    residual!,
    jacobian!,
    mass_matrix,

    NullOperator,
    OperatorSum,
    active_operators,

    AbstractDiffusionCoefficients,
    get_D,
    ConstantDiffusion,
    ArrheniusDiffusion,
    CallableDiffusion,

    LinearDiffusionOperator,
    WeakDirichletBoundaryOperator,
    DirichletBoundaryOperator,
    PenaltyMethod,
    MassMatrixMethod,
    CallbackMethod,
    EliminatedMethod,
    build_solver_callback,
    ToyReactionOperator,
    SimpleTrappingReactionOperator,
    ConstraintOperator,

    UnsplitFormulation,
    IMEXFormulation,
    IMEXReactionFormulation,
    LieSplit,
    StrangSplit,
    SplitFormulation,
    ResidualFormulation,

    build_rd_model,
    build_trapping_model,
    trapping_variable_layout,

    SplitProblem,
    SplitSolution,

    build_problem,
    solve_problem,

    load_config,
    run_simulation,

    SimulationResult,
    wrap_result,
    variable_timeseries,
    variable_snapshot,
    integrated_variable,
    surface_diffusive_fluxes,
    build_summary_dataframe,
    write_field_output_hdf5,
    write_summary_csv,
    load_ic_from_hdf5,
    print_run_banner,
    solver_stats_dict,
    check_mass_conservation,

    AbstractTemperatureProvider,
    ConstantTemperature,
    LinearRampTemperature,
    FunctionTemperature,
    temperature_at,

    HotgatesTrappingAdaptor,
    HotgatesReactionOperator,
    build_hotgates_variable_layout,
    build_hotgates_trapping_model,
    build_palioxis_trapping_model,
    build_equilibrium_ic,
    build_ic_from_total_hydrogen,

    FakeHotgatesModel,
    hotgates_rates!,

    plot_tds_flux,
    plot_spatial_snapshot,
    plot_spatial_evolution,
    record_spatial_video,
    record_spatial_animation
end
