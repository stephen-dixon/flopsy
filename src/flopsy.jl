module Flopsy

using LinearAlgebra
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

    LinearDiffusionOperator,
    ToyReactionOperator,
    SimpleTrappingReactionOperator,
    ConstraintOperator,

    UnsplitFormulation,
    IMEXFormulation,
    LieSplit,
    StrangSplit,
    SplitFormulation,
    ResidualFormulation,

    build_rd_model,
    build_trapping_model,
    trapping_variable_layout,

    build_problem,
    solve_problem,

    load_config,
    run_simulation,

    HotgatesReactionOperator

    SimulationResult,
    wrap_result,
    variable_timeseries,
    variable_snapshot,
    integrated_variable,
    build_summary_dataframe,
    write_field_output_hdf5,
    write_summary_csv,
    print_run_banner,
    solver_stats_dict
end
