"""
    validate(cfg::ProblemConfig) -> ProblemConfig

Validate a parsed configuration and throw user-facing errors before solver
construction.
"""
function validate(cfg::ProblemConfig)
    cfg.mesh.kind == :uniform_1d ||
        throw(ArgumentError("Only mesh.kind = \"uniform_1d\" is currently supported"))
    cfg.mesh.nx >= 2 || throw(ArgumentError("mesh.nx must be at least 2"))
    cfg.mesh.xmax > cfg.mesh.xmin ||
        throw(ArgumentError("mesh.xmax must be greater than mesh.xmin"))

    cfg.problem_type in (:diffusion_1d, :trapping_1d, :hotgates_trapping) ||
        throw(ArgumentError("Unsupported problem_type $(cfg.problem_type). Supported types: diffusion_1d, trapping_1d, hotgates_trapping"))

    cfg.solver.abstol > 0 || throw(ArgumentError("solver.abstol must be positive"))
    cfg.solver.reltol > 0 || throw(ArgumentError("solver.reltol must be positive"))

    t0 = getproperty(cfg.parameters, :t0)
    tend = getproperty(cfg.parameters, :tend)
    tend > t0 || throw(ArgumentError("parameters.tend must be greater than parameters.t0"))

    if cfg.solver.formulation == :split
        cfg.solver.dt === nothing &&
            throw(ArgumentError("solver.dt is required for split formulation"))
        cfg.solver.dt > 0 || throw(ArgumentError("solver.dt must be positive"))
    end

    if cfg.solver.saveat !== nothing
        issorted(cfg.solver.saveat) ||
            throw(ArgumentError("solver.saveat must be sorted in ascending order"))
        for t in cfg.solver.saveat
            t0 <= t <= tend ||
                throw(ArgumentError("solver.saveat time $t lies outside [$t0, $tend]"))
        end
    end

    for bc in cfg.boundary_conditions
        bc.side in (:left, :right) ||
            throw(ArgumentError("boundary condition side must be :left or :right"))
        bc.kind == :dirichlet ||
            throw(ArgumentError("Only dirichlet boundary conditions are currently supported"))
        bc.method in (:weak, :penalty, :mass_matrix, :callback, :eliminated) ||
            throw(ArgumentError("Unsupported boundary condition method $(bc.method)"))
    end

    return cfg
end
