"""
    SplitProblem

Internal problem type returned by `build_problem` for `SplitFormulation`.
Holds the model, initial condition, time span, and solver configuration so that
`solve_problem` can execute the splitting loop.
"""
struct SplitProblem
    model::SystemModel
    u0::Vector{Float64}
    tspan::Tuple{Float64,Float64}
    solver_config::SolverConfig
end

"""
    SplitSolution

Minimal solution object returned by the operator-splitting integrator.
Provides `.t` and `.u` fields matching the interface expected by Flopsy output
helpers (`variable_timeseries`, `build_summary_dataframe`, etc.).
`retcode` is `:Success` when the integration completed without error.
"""
struct SplitSolution
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
    retcode::Symbol
end

# Make SplitSolution behave like a SciML solution for the subset of fields
# used by Flopsy output helpers.
Base.length(sol::SplitSolution) = length(sol.t)


"""
    build_problem(model, u0, tspan, ::SplitFormulation, solver_config) -> SplitProblem

Return a `SplitProblem` that encapsulates everything needed for the splitting
loop.  The actual integration is performed by `solve_problem(::SplitProblem, ...)`.
"""
function build_problem(model::SystemModel, u0, tspan, ::SplitFormulation, solver_config::SolverConfig)
    return SplitProblem(model, collect(Float64, u0),
                        (Float64(tspan[1]), Float64(tspan[2])),
                        solver_config)
end


# ---------------------------------------------------------------------------
# Internal splitting loop
# ---------------------------------------------------------------------------

function _solve_split(prob::SplitProblem, solver_config::SolverConfig, ::LieSplit)
    t0, t_end = prob.tspan
    Δt        = Float64(solver_config.dt)
    model     = prob.model
    ctx       = model.context

    react_ops = Tuple(filter(!isnothing, [model.operators.reaction]))
    diff_ops  = Tuple(filter(!isnothing, [model.operators.diffusion,
                                          model.operators.boundary]))

    react_op = isempty(react_ops) ? NullOperator() : OperatorSum(react_ops)
    diff_op  = isempty(diff_ops)  ? NullOperator() : OperatorSum(diff_ops)

    saveat_set = _build_saveat_set(solver_config, t0, t_end)

    t = t0
    u = copy(prob.u0)

    saved_t = Float64[]
    saved_u = Vector{Float64}[]

    _record!(saved_t, saved_u, t, u, saveat_set)

    while t < t_end - sqrt(eps(t_end))
        dt = min(Δt, t_end - t)

        # Lie: reaction sub-step → diffusion sub-step
        u = _solve_substep(react_op, ctx, u, t, dt, solver_config)
        u = _solve_substep(diff_op,  ctx, u, t, dt, solver_config)

        t += dt
        _record!(saved_t, saved_u, t, u, saveat_set)
    end

    return SplitSolution(saved_t, saved_u, :Success)
end

function _solve_split(prob::SplitProblem, solver_config::SolverConfig, ::StrangSplit)
    t0, t_end = prob.tspan
    Δt        = Float64(solver_config.dt)
    model     = prob.model
    ctx       = model.context

    react_ops = Tuple(filter(!isnothing, [model.operators.reaction]))
    diff_ops  = Tuple(filter(!isnothing, [model.operators.diffusion,
                                          model.operators.boundary]))

    react_op = isempty(react_ops) ? NullOperator() : OperatorSum(react_ops)
    diff_op  = isempty(diff_ops)  ? NullOperator() : OperatorSum(diff_ops)

    saveat_set = _build_saveat_set(solver_config, t0, t_end)

    t = t0
    u = copy(prob.u0)

    saved_t = Float64[]
    saved_u = Vector{Float64}[]

    _record!(saved_t, saved_u, t, u, saveat_set)

    while t < t_end - sqrt(eps(t_end))
        dt = min(Δt, t_end - t)

        # Strang: half reaction → full diffusion → half reaction
        u = _solve_substep(react_op, ctx, u, t,          dt / 2, solver_config)
        u = _solve_substep(diff_op,  ctx, u, t,          dt,     solver_config)
        u = _solve_substep(react_op, ctx, u, t + dt / 2, dt / 2, solver_config)

        t += dt
        _record!(saved_t, saved_u, t, u, saveat_set)
    end

    return SplitSolution(saved_t, saved_u, :Success)
end


# ---------------------------------------------------------------------------
# Sub-step helpers
# ---------------------------------------------------------------------------

"""Solve one operator sub-problem over [t, t+dt] starting from u."""
function _solve_substep(op, ctx::SystemContext, u::Vector{Float64},
                        t::Float64, dt::Float64,
                        solver_config::SolverConfig)
    op isa NullOperator && return u

    function f!(du, uv, p, tv)
        rhs!(du, op, uv, ctx, tv)
        return nothing
    end

    ode_f = ODEFunction(f!)

    sub_prob = ODEProblem(ode_f, copy(u), (t, t + dt))

    sub_sol = SciMLBase.solve(
        sub_prob,
        solver_config.algorithm;
        abstol         = solver_config.abstol,
        reltol         = solver_config.reltol,
        save_everystep = false,
    )

    return copy(sub_sol.u[end])
end


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

"""Return a Set of save times within [t0, t_end] from solver_config.saveat."""
function _build_saveat_set(solver_config::SolverConfig, t0::Float64, t_end::Float64)
    saveat = solver_config.saveat
    saveat === nothing && return Set{Float64}()
    return Set{Float64}(filter(s -> t0 <= s <= t_end, saveat))
end

"""Record (t, u) if t matches a save time or no explicit saveat is given."""
function _record!(saved_t, saved_u, t::Float64, u::Vector{Float64},
                  saveat_set::Set{Float64})
    should_save = if isempty(saveat_set)
        true
    else
        any(s -> abs(t - s) <= sqrt(eps(max(abs(t), abs(s), 1.0))), saveat_set)
    end

    if should_save
        # Avoid duplicate entries (e.g. t0 recorded twice at loop boundary).
        if isempty(saved_t) || saved_t[end] != t
            push!(saved_t, t)
            push!(saved_u, copy(u))
        end
    end
end
