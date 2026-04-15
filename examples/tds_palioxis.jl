# =============================================================================
# TDS example: hydrogen implantation and thermal desorption spectroscopy
# using a real Palioxis MultipleDefectModel (via PalioxisExt).
#
# Physics
# -------
#   - 1 mobile (diffusing) hydrogen species
#   - 1 defect type (vacancy), uniform concentration 0.001 everywhere
#   - 6 non-equilibrium trapped species (occupancy levels 1–6 in the vacancy)
#
# Simulation stages
# -----------------
#   1. Loading   (0 → 24 h)    : T = 400 K, left BC = 1.0, right BC = 0 (vacuum)
#   2. Rest      (24 → 26 h)   : T = 300 K, both BCs = 0
#   3. TDS ramp  (26 h → end)  : T ramps 300 K → 1200 K at 10 K/min, both BCs = 0
#
# The diffusion coefficient D(T) is evaluated from the Palioxis model at every
# time step — fully temperature-dependent, consistent with the model parameters.
#
# Outputs
# -------
#   tds_summary.csv  — time, temperature, integrated inventories, surface fluxes
#   tds_fields.h5    — full (nt, nx) field output for all variables
# =============================================================================

using Flopsy
using Palioxis

using LinearAlgebra
using SparseArrays
using ADTypes
using TOML
using SciMLBase
using OrdinaryDiffEq

# ---------------------------------------------------------------------------
# 1.  Palioxis model
# ---------------------------------------------------------------------------

const PALIOXIS_ROOT = get(ENV, "PALIOXIS_ROOT", "")
isempty(PALIOXIS_ROOT) || Palioxis.init(PALIOXIS_ROOT)
const XML_FILE = get(ENV, "PALIOXIS_XML", "traps.xml")

pal = MultipleDefectModel(XML_FILE)
report(pal)

@assert get_n_gas(pal) == 1 "expected 1 mobile species"
@assert get_n_trap_types(pal) == 1 "expected 1 defect/trap type"
@assert get_n_ne_traps(pal) == 6 "expected 6 occupancy-level DOFs"

println("Mobile species : ", gas_names(pal))
println("Trap species   : ", trap_names(pal))
println("Defect names   : ", defect_names(pal))

# ---------------------------------------------------------------------------
# 2.  Spatial mesh
# ---------------------------------------------------------------------------

const L = 1.0e-3   # sample thickness [m]
const NX = 200

mesh = Mesh1D(L, NX)

# ---------------------------------------------------------------------------
# 3.  Defect profile (uniform)
# ---------------------------------------------------------------------------

const RHO_DEFECT = 0.001
defects = fill(RHO_DEFECT, get_n_trap_types(pal), NX)   # (n_traps, nx)

# ---------------------------------------------------------------------------
# 4.  Time structure
# ---------------------------------------------------------------------------

const t_load = 24.0 * 3600.0
const t_rest_end = t_load + 2.0 * 3600.0

const T_RAMP_START = 300.0
const T_RAMP_END = 1200.0
const RAMP_RATE = 10.0 / 60.0           # K/s  (10 K/min)
const t_tds_end = t_rest_end + (T_RAMP_END - T_RAMP_START) / RAMP_RATE

const tspan = (0.0, t_tds_end)

# ---------------------------------------------------------------------------
# 5.  Temperature provider  (piecewise: load / rest / ramp)
# ---------------------------------------------------------------------------

function T_profile(t::Real)::Float64
    if t < t_load
        return 400.0
    elseif t < t_rest_end
        return 300.0
    else
        return min(T_RAMP_START + RAMP_RATE * (t - t_rest_end), T_RAMP_END)
    end
end

temperature = FunctionTemperature(T_profile)

# ---------------------------------------------------------------------------
# 6.  Build model
#     D(T) comes from Palioxis at every step — no pre-evaluation needed.
#     Boundary conditions are passed as callables: the extension constructs
#     DirichletBoundaryOperator internally, matched to the same D(T).
# ---------------------------------------------------------------------------

const C_SURFACE_LOAD = 1.0e-10

model = build_palioxis_trapping_model(
    palioxis_model = pal,
    mesh = mesh,
    defects = defects,
    temperature = temperature,
    left_bc = t -> t < t_load ? C_SURFACE_LOAD : 0.0,
    right_bc = t -> 0.0
)

println("\nVariable layout:")
for (i, name) in enumerate(variable_names(model.layout))
    println("  [$i] $name")
end

# ---------------------------------------------------------------------------
# 7.  Initial condition
#
#     Option A — all zeros, let the loading BC drive hydrogen in.
#                Simple; loading phase transient is unphysical for ~first dx/D.
#
#     Option B — equilibrium IC from a user-supplied total H profile.
#                Use when you have a known implantation damage profile.
# ---------------------------------------------------------------------------

# Option A: start from a hydrogen-free sample
u0 = zeros(Float64, nvariables(model.layout) * NX)

# Set the left boundary node to match the loading BC at t=0 to avoid a
# step discontinuity.
state_view(u0, model.layout, NX)[1, 1] = C_SURFACE_LOAD

# Option B (uncomment to use):
# # Gaussian implantation profile peaked at 50 µm depth, σ = 20 µm
# x  = model.context.mesh.x
# C0 = 0.5 .* exp.(.-((x .- 50e-6).^2) ./ (2 * (20e-6)^2))
# u0 = build_ic_from_total_hydrogen(pal, model, C0, 400.0)

# ---------------------------------------------------------------------------
# 8.  Solver configuration
# ---------------------------------------------------------------------------

saveat = sort(unique(vcat(
    collect(range(0.0, t_load, step = 1800.0)),
    collect(range(t_load, t_rest_end, step = 300.0)),
    collect(range(t_rest_end, t_tds_end, step = 10.0)),
    [t_tds_end]
)))

solver_config = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm = Rodas5(autodiff = AutoFiniteDiff()),
    abstol = 1e-10,
    reltol = 1e-8,
    saveat = saveat
)

print_run_banner(Dict{String, Any}(), solver_config, model)

# ---------------------------------------------------------------------------
# 9.  Solve
# ---------------------------------------------------------------------------

t_wall = @elapsed sol = solve_problem(model, u0, tspan, solver_config)
println("Solve complete: retcode = $(sol.retcode)  wall time = $(round(t_wall; digits=2)) s")

# ---------------------------------------------------------------------------
# 10.  Post-processing
# ---------------------------------------------------------------------------

# Compute temperature and total inventory at every saved time for the CSV.
t_saved = sol.t
T_saved = T_profile.(t_saved)

tmp_result = wrap_result(model, sol, nothing)
all_names = variable_names(model.layout)
total_inv = sum(integrated_variable(tmp_result, name) for name in all_names)

result = wrap_result(
    model, sol, nothing;
    summaries = Dict{Symbol, Any}(
        :extra_timeseries => Dict(
        "temperature_K" => T_saved,
        "total_inventory" => total_inv
    ),
    ),
    metadata = Dict{String, Any}(
        "model_type" => "tds_palioxis",
        "xml_file" => XML_FILE,
        "L_m" => L,
        "nx" => NX,
        "rho_defect" => RHO_DEFECT,
        "T_loading_K" => 400.0,
        "T_rest_K" => 300.0,
        "ramp_rate_K_per_s" => RAMP_RATE,
        "retcode" => string(sol.retcode),
        "walltime_s" => t_wall
    )
)

# Summary CSV columns:
#   time, temperature_K, total_inventory,
#   integral_<var>  × (1 + 6),
#   left_flux_<mobile_var>, right_flux_<mobile_var>
write_summary_csv(result, "tds_summary.csv")
println("Written: tds_summary.csv")

write_field_output_hdf5(result, "tds_fields.h5")
println("Written: tds_fields.h5")

# ---------------------------------------------------------------------------
# 11.  Quick diagnostics: peak desorption temperature
# ---------------------------------------------------------------------------

fluxes = surface_diffusive_fluxes(result)
mobile_idx = first(variables_in_group(model.layout, :mobile))
mobile_sym = variable_names(model.layout)[mobile_idx]

i_ramp_start = findfirst(>=(t_rest_end), t_saved)
right_flux_ramp = fluxes[mobile_sym].right[i_ramp_start:end]
T_ramp = T_saved[i_ramp_start:end]

i_peak = argmax(right_flux_ramp)
println("\nPeak right-surface desorption flux during TDS ramp:")
println("  flux = $(right_flux_ramp[i_peak])")
println("  T    = $(T_ramp[i_peak]) K")
println("  t    = $(round((t_saved[i_ramp_start-1+i_peak] - t_rest_end)/60; digits=1)) min into ramp")
