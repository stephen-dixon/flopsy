"""
Formulation comparison — same 1-D reaction-diffusion trapping problem solved with
four different Flopsy formulations:

  1. UnsplitFormulation   — monolithic stiff ODE with analytic sparse Jacobian
  2. IMEXFormulation      — IMEX split: diffusion implicit, reaction explicit
  3. SplitFormulation(LieSplit)    — operator splitting (first-order Lie)
  4. SplitFormulation(StrangSplit) — operator splitting (second-order Strang)
  5. ResidualFormulation  — DAE with mass matrix (quasi-static trapping constraint)

Run from the flopsy root with:
    julia --project=. examples/formulation_comparison.jl
"""

using Flopsy
using OrdinaryDiffEq
using ADTypes: AutoFiniteDiff
using LinearAlgebra: norm

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------

mesh = Mesh1D(1e-3, 40)          # 1 mm domain, 40 nodes
nx = length(mesh.x)
t_end = 200.0
tspan = (0.0, t_end)
n_save = 101
saveat = range(0.0, t_end; length = n_save)

k_trap = 5e-3
k_detrap = 0.05
D_mobile = 1e-7                     # m²/s

model = build_trapping_model(
    mesh = mesh,
    k_trap = k_trap,
    k_detrap = k_detrap,
    diffusion_coefficient = D_mobile
)

# Gaussian initial mobile profile; traps empty
x = mesh.x
c0 = 0.05 .* exp.(-(x .- x[end]/2) .^ 2 ./ (2*(x[end]/8)^2))
th0 = zeros(nx)

# Pack into the flat state vector (nvars=2, column-major: [c, theta] per node)
u0 = zeros(2 * nx)
for ix in 1:nx
    u0[(ix - 1) * 2 + 1] = c0[ix]
    u0[(ix - 1) * 2 + 2] = th0[ix]
end

println("Model: nx=$(nx)  t_end=$(t_end)  k_trap=$(k_trap)  k_detrap=$(k_detrap)  D=$(D_mobile)")
println()

# ---------------------------------------------------------------------------
# Helper: extract mobile and trap profiles from a SimulationResult
# ---------------------------------------------------------------------------

function extract_profiles(result)
    c = variable_timeseries(result, :c)
    th = variable_timeseries(result, :theta)
    return c, th
end

# ---------------------------------------------------------------------------
# 1. UnsplitFormulation (reference)
# ---------------------------------------------------------------------------

println("=== 1. UnsplitFormulation ===")
config_unsplit = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm = Rodas5P(),
    abstol = 1e-9,
    reltol = 1e-7,
    saveat = collect(saveat),
    show_progress = false,
    show_solver_stats = false
)

t_unsplit = @elapsed begin
    sol_unsplit = solve_problem(model, u0, tspan, config_unsplit)
    result_unsplit = wrap_result(model, sol_unsplit, config_unsplit)
end
c_ref, th_ref = extract_profiles(result_unsplit)
println("  Saved $(length(result_unsplit.solution.t)) time points  |  wall time: $(round(t_unsplit; digits=2)) s")
println("  Final total-H: $(round(sum(c_ref[end,:] .+ th_ref[end,:]) * mesh.dx; sigdigits=5))")
println()

# ---------------------------------------------------------------------------
# 2. IMEXFormulation
# ---------------------------------------------------------------------------

println("=== 2. IMEXFormulation (KenCarp4) ===")
config_imex = SolverConfig(
    formulation = IMEXFormulation(),
    algorithm = KenCarp4(),
    abstol = 1e-9,
    reltol = 1e-7,
    saveat = collect(saveat),
    show_progress = false,
    show_solver_stats = false
)

t_imex = @elapsed begin
    sol_imex = solve_problem(model, u0, tspan, config_imex)
    result_imex = wrap_result(model, sol_imex, config_imex)
end
c_imex, th_imex = extract_profiles(result_imex)
err_imex = norm(c_imex .- c_ref) / (norm(c_ref) + 1e-30)
println("  Saved $(length(result_imex.solution.t)) time points  |  wall time: $(round(t_imex; digits=2)) s")
println("  Relative L2 error vs. Unsplit (mobile): $(round(err_imex; sigdigits=4))")
println()

# ---------------------------------------------------------------------------
# 3. SplitFormulation — LieSplit
# ---------------------------------------------------------------------------

println("=== 3. SplitFormulation(LieSplit) ===")
# Sub-steps have no analytic Jacobian; ForwardDiff Dual numbers conflict with
# pre-allocated Float64 scratch buffers — use AutoFiniteDiff instead.
config_lie = SolverConfig(
    formulation = SplitFormulation(LieSplit()),
    algorithm = Rodas5(autodiff = AutoFiniteDiff()),
    abstol = 1e-9,
    reltol = 1e-7,
    saveat = collect(saveat),
    dt = 2.0,             # macro-step size (seconds)
    show_progress = false,
    show_solver_stats = false
)

t_lie = @elapsed begin
    sol_lie = solve_problem(model, u0, tspan, config_lie)
    result_lie = wrap_result(model, sol_lie, config_lie)
end
c_lie, th_lie = extract_profiles(result_lie)

# Compare at shared times only
shared_lie = [findfirst(t -> abs(t - s) < 0.5 * config_lie.dt, sol_lie.t)
              for s in result_unsplit.solution.t]
shared_lie = filter(!isnothing, shared_lie)
if !isempty(shared_lie)
    err_lie = norm(c_lie[shared_lie, :] .- c_ref[1:length(shared_lie), :]) /
              (norm(c_ref[1:length(shared_lie), :]) + 1e-30)
    println("  Saved $(length(sol_lie.t)) time points  |  wall time: $(round(t_lie; digits=2)) s")
    println("  Relative L2 error vs. Unsplit (mobile): $(round(err_lie; sigdigits=4))")
else
    println("  Saved $(length(sol_lie.t)) time points  |  wall time: $(round(t_lie; digits=2)) s")
end
println()

# ---------------------------------------------------------------------------
# 4. SplitFormulation — StrangSplit
# ---------------------------------------------------------------------------

println("=== 4. SplitFormulation(StrangSplit) ===")
config_strang = SolverConfig(
    formulation = SplitFormulation(StrangSplit()),
    algorithm = Rodas5(autodiff = AutoFiniteDiff()),
    abstol = 1e-9,
    reltol = 1e-7,
    saveat = collect(saveat),
    dt = 2.0,
    show_progress = false,
    show_solver_stats = false
)

t_strang = @elapsed begin
    sol_strang = solve_problem(model, u0, tspan, config_strang)
    result_strang = wrap_result(model, sol_strang, config_strang)
end
c_strang, th_strang = extract_profiles(result_strang)
println("  Saved $(length(sol_strang.t)) time points  |  wall time: $(round(t_strang; digits=2)) s")
println()

# ---------------------------------------------------------------------------
# 5. ResidualFormulation (DAE / quasi-static trapping)
# ---------------------------------------------------------------------------

println("=== 5. ResidualFormulation (Rodas5P, mass-matrix DAE) ===")
config_dae = SolverConfig(
    formulation = ResidualFormulation(),
    algorithm = Rodas5P(),
    abstol = 1e-9,
    reltol = 1e-7,
    saveat = collect(saveat),
    show_progress = false,
    show_solver_stats = false
)

t_dae = @elapsed begin
    sol_dae = solve_problem(model, u0, tspan, config_dae)
    result_dae = wrap_result(model, sol_dae, config_dae)
end
c_dae, th_dae = extract_profiles(result_dae)
err_dae = norm(c_dae .- c_ref) / (norm(c_ref) + 1e-30)
println("  Saved $(length(result_dae.solution.t)) time points  |  wall time: $(round(t_dae; digits=2)) s")
println("  Relative L2 error vs. Unsplit (mobile): $(round(err_dae; sigdigits=4))")
println()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

println("=== Summary ===")
println("  UnsplitFormulation :  $(round(t_unsplit; digits=2)) s  (reference)")
println("  IMEXFormulation    :  $(round(t_imex;   digits=2)) s  |  L2 err ≈ $(round(err_imex; sigdigits=3))")
println("  LieSplit           :  $(round(t_lie;    digits=2)) s")
println("  StrangSplit        :  $(round(t_strang; digits=2)) s")
println("  ResidualFormulation:  $(round(t_dae;    digits=2)) s  |  L2 err ≈ $(round(err_dae;  sigdigits=3))")
println()
println("Notes:")
println("  - Unsplit and IMEX use adaptive stepping; Split uses a fixed macro-step")
println("    of dt=$(config_lie.dt) s with sub-step tolerance 1e-9/1e-7.")
println("  - The ResidualFormulation (DAE) enforces quasi-static trap equilibrium")
println("    (mass-matrix=0 for :trap variables), which is a *different physical model*")
println("    from full dynamic trapping. Error vs Unsplit reflects the quasi-static")
println("    approximation, not solver error — appropriate only when trapping is fast.")
println("  - SplitFormulation sub-steps use AutoFiniteDiff because pre-allocated")
println("    Float64 scratch buffers are incompatible with ForwardDiff Dual numbers.")
