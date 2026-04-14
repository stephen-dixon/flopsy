# =============================================================================
# TDS from equilibrium initial condition
#
# Physics
# -------
#   - 1 mobile (diffusing) hydrogen species
#   - 1 defect type (vacancy), uniform concentration 0.001
#   - 6 non-equilibrium trapped species (occupancy levels 1-6)
#
# Approach
# --------
#   Start from an equilibrium IC corresponding to a uniform total H
#   concentration of 1e-10 at 300 K, computed by partitioning into
#   mobile + trapped via Palioxis.calculate_steady_state!.
#   No loading or rest phase — only the TDS ramp (300 → 1200 K, 10 K/min).
#
# Outputs
# -------
#   out_eq/tds_eq_summary.csv
#   out_eq/tds_eq_fields.h5
#   out_eq/fig_tds_flux.png
#   out_eq/fig_spatial_initial.png
#   out_eq/fig_spatial_final.png
#   out_eq/spatial_animation.mp4
# =============================================================================

using Flopsy
using Palioxis

using LinearAlgebra
using SparseArrays
using ADTypes
using SciMLBase
using OrdinaryDiffEq
using CairoMakie
using Printf

mkpath("out_eq")

# ---------------------------------------------------------------------------
# 1.  Palioxis model
# ---------------------------------------------------------------------------

Palioxis.init(get(ENV, "PALIOXIS_ROOT", "/Users/sdixon/src/palioxis-tds/Palioxis"))
const XML_FILE = get(ENV, "PALIOXIS_XML", "traps.xml")

pal = MultipleDefectModel(XML_FILE)
report(pal)

@assert get_n_gas(pal)        == 1  "expected 1 mobile species"
@assert get_n_trap_types(pal) == 1  "expected 1 defect/trap type"
@assert get_n_ne_traps(pal)   == 6  "expected 6 occupancy-level DOFs"

println("Mobile species : ", gas_names(pal))
println("Trap species   : ", trap_names(pal))

# ---------------------------------------------------------------------------
# 2.  Spatial mesh
# ---------------------------------------------------------------------------

const L  = 1.0e-3   # sample thickness [m]
const NX = 200

mesh = Mesh1D(L, NX)

# ---------------------------------------------------------------------------
# 3.  Defect profile
# ---------------------------------------------------------------------------

const RHO_DEFECT = 0.001
defects = fill(RHO_DEFECT, get_n_trap_types(pal), NX)

# ---------------------------------------------------------------------------
# 4.  TDS ramp: 300 K → 1200 K at 10 K/min, vacuum BCs throughout
# ---------------------------------------------------------------------------

const T_START   = 300.0
const T_END     = 1200.0
const RAMP_RATE = 10.0 / 60.0          # K/s
const t_ramp    = (T_END - T_START) / RAMP_RATE
const tspan     = (0.0, t_ramp)

T_profile(t::Real)::Float64 = min(T_START + RAMP_RATE * t, T_END)
temperature = FunctionTemperature(T_profile)

# ---------------------------------------------------------------------------
# 5.  Build model with vacuum (zero-concentration) BCs on both surfaces
# ---------------------------------------------------------------------------

model = build_palioxis_trapping_model(
    palioxis_model = pal,
    mesh           = mesh,
    defects        = defects,
    temperature    = temperature,
    left_bc        = t -> 0.0,
    right_bc       = t -> 0.0,
)

println("\nVariable layout:")
for (i, name) in enumerate(variable_names(model.layout))
    println("  [$i] $name")
end

# ---------------------------------------------------------------------------
# 6.  Equilibrium IC from uniform total H = 1e-10 at 300 K
# ---------------------------------------------------------------------------

const C_TOT = 1e-10
total_H     = fill(C_TOT, NX)

println("\nBuilding equilibrium IC from uniform total H = $C_TOT at $(T_START) K ...")
u0 = build_ic_from_total_hydrogen(pal, model, total_H, T_START)

# Ensure boundary nodes match the vacuum BC to avoid step discontinuity.
nvars = nvariables(model.layout)
U0    = reshape(u0, nvars, NX)
mobile_idx = first(variables_in_group(model.layout, :mobile))
U0[mobile_idx, 1]  = 0.0
U0[mobile_idx, NX] = 0.0

println("IC built.")

# Print initial inventories for sanity check.
println("Initial total mobile inventory (node sum) : ",
        sum(U0[mobile_idx, :]) * mesh.dx)

# ---------------------------------------------------------------------------
# 7.  Save-points: dense during ramp
# ---------------------------------------------------------------------------

saveat = collect(range(0.0, t_ramp; step = 10.0))
push!(saveat, t_ramp)
unique!(sort!(saveat))

# ---------------------------------------------------------------------------
# 8.  Solver configuration (UnsplitFormulation with Rodas5)
# ---------------------------------------------------------------------------

solver_config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm   = Rodas5(autodiff = AutoFiniteDiff()),
    abstol      = 1e-10,
    reltol      = 1e-8,
    saveat      = saveat,
)

print_run_banner(Dict{String,Any}(), solver_config, model)

# ---------------------------------------------------------------------------
# 9.  Solve
# ---------------------------------------------------------------------------

t_wall = @elapsed sol = solve_problem(model, u0, tspan, solver_config)
println("Solve complete: retcode = $(sol.retcode)  wall time = $(round(t_wall; digits=2)) s")

# ---------------------------------------------------------------------------
# 10.  Wrap result with temperature timeseries
# ---------------------------------------------------------------------------

t_saved = sol.t
T_saved = T_profile.(t_saved)

all_names = variable_names(model.layout)

result = wrap_result(
    model, sol, nothing;
    summaries = Dict{Symbol,Any}(
        :extra_timeseries => Dict(
            "temperature_K" => T_saved,
        ),
    ),
    metadata = Dict{String,Any}(
        "model_type"        => "tds_from_equilibrium",
        "xml_file"          => XML_FILE,
        "L_m"               => L,
        "nx"                => NX,
        "rho_defect"        => RHO_DEFECT,
        "C_tot_IC"          => C_TOT,
        "T_start_K"         => T_START,
        "T_end_K"           => T_END,
        "ramp_rate_K_per_s" => RAMP_RATE,
        "retcode"           => string(sol.retcode),
        "walltime_s"        => t_wall,
    ),
)

write_summary_csv(result, "out_eq/tds_eq_summary.csv")
write_field_output_hdf5(result, "out_eq/tds_eq_fields.h5")
println("Written: out_eq/tds_eq_summary.csv, out_eq/tds_eq_fields.h5")

# ---------------------------------------------------------------------------
# 11.  Peak desorption temperature
# ---------------------------------------------------------------------------

fluxes     = surface_diffusive_fluxes(result)
mobile_sym = variable_names(model.layout)[mobile_idx]
right_flux = fluxes[mobile_sym].right

i_peak = argmax(right_flux)
println("\nPeak right-surface desorption flux:")
println("  flux = $(right_flux[i_peak])")
println("  T    = $(T_saved[i_peak]) K")

# ---------------------------------------------------------------------------
# 12.  Grouping helpers for plots
# ---------------------------------------------------------------------------

# Group mobile and sum all trap occupancy levels into "Total trapped".
group_fn = (names, data) -> begin
    out = Dict{String,Vector{Float64}}()
    mob_i  = [i for (i,n) in enumerate(names) if startswith(String(n), "mobile")]
    trap_i = [i for (i,n) in enumerate(names) if startswith(String(n), "trap")]
    for i in mob_i
        out[String(names[i])] = vec(copy(data[i, :]))
    end
    if !isempty(trap_i)
        tot = zeros(size(data, 2))
        for i in trap_i; tot .+= data[i, :]; end
        out["Total trapped"] = tot
    end
    out
end

# All individual variables (for the animation — 1 mobile + 6 trap levels = 7 lines).
group_all = (names, data) -> Dict(String(n) => vec(copy(data[i,:])) for (i,n) in enumerate(names))

# ---------------------------------------------------------------------------
# 13.  Spatial plots
# ---------------------------------------------------------------------------

# TDS flux vs temperature.
fig_flux = plot_tds_flux(result;
    surface = :right,
    title   = "TDS desorption spectrum (from equilibrium IC)",
)
save("out_eq/fig_tds_flux.png", fig_flux)

# Initial spatial profile (all 7 fields on one log-log axis).
fig_init = plot_spatial_snapshot(result;
    time_index = 1,
    group_fn   = group_all,
    all_on_one = true,
    xscale     = :log10,
    yscale     = :log10,
    xmin       = 1e-20,
    ymin       = 1e-20,
    title      = @sprintf("Initial spatial profile  (T = %.0f K)", T_START),
)
save("out_eq/fig_spatial_initial.png", fig_init)

# Final spatial profile.
fig_final = plot_spatial_snapshot(result;
    time_index = :last,
    group_fn   = group_all,
    all_on_one = true,
    xscale     = :log10,
    yscale     = :log10,
    xmin       = 1e-20,
    ymin       = 1e-20,
    title      = @sprintf("Final spatial profile  (T = %.0f K)", T_saved[end]),
)
save("out_eq/fig_spatial_final.png", fig_final)

println("Written: out_eq/fig_tds_flux.png, out_eq/fig_spatial_initial.png, out_eq/fig_spatial_final.png")

# ---------------------------------------------------------------------------
# 14.  Animation — all 7 fields (1 mobile + 6 trap levels) on one log-log axis
# ---------------------------------------------------------------------------

println("Rendering animation ($(length(t_saved)) frames) ...")
anim_path = record_spatial_animation(
    result,
    "out_eq/spatial_animation.mp4";
    group_fn = group_all,
    fps      = 20,
    xscale   = :log10,
    yscale   = :log10,
    xmin     = 1e-20,
    ymin     = 1e-20,
    colormap = :tab10,
    title_fn = t -> begin
        T_t = T_profile(t)
        @sprintf("t = %.0f s   T = %.1f K", t, T_t)
    end,
)
println("Written: $anim_path")

println("\nDone.")
