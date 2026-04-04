module Flopsy

using LinearAlgebra

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

include("runner.jl")
include("output.jl")

export
    # core types
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

    # layout/state helpers
    nvariables,
    variable_names,
    variables_in_group,
    variables_with_tag,
    has_group,
    state_view,
    node_view,
    variable_view,
    group_view,

    # operator capability hooks
    supports_rhs,
    supports_implicit_rhs,
    supports_step,
    supports_residual,
    supports_mass_matrix,
    supports_jacobian,

    # operator API
    rhs!,
    implicit_rhs!,
    step!,
    residual!,
    jacobian!,
    mass_matrix,

    # operator utilities
    NullOperator,
    OperatorSum,
    active_operators,

    # concrete starter operators
    LinearDiffusionOperator,
    ConstraintOperator,

    # formulations
    UnsplitFormulation,
    IMEXFormulation,
    LieSplit,
    StrangSplit,
    SplitFormulation,
    ResidualFormulation,

    # model builders
    build_rd_model,
    build_trapping_model,

    # problem / solver
    build_problem,
    solve_problem,

    # config / runner
    load_config,
    run_simulation,

    # adaptor-facing helper
    HotgatesReactionOperator

end
