"""
Package extension that wires the real `Palioxis.MultipleDefectModel` into
Flopsy's Hotgates adapter.  Activated automatically when both `Flopsy` and
`Palioxis` are loaded in the same Julia session.
"""
module PalioxisExt

using Flopsy
using Palioxis
import CSV
import DataFrames

function register_flopsy_plugin!(registry::Flopsy.SyntaxRegistry)
    Flopsy.register_syntax!(registry,
        Flopsy.SyntaxSpec(
            :backend,
            :palioxis_trapping,
            [
                Flopsy.ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
                Flopsy.ParameterSpec(
                    :xml_file, true, nothing, "Palioxis XML model file"; kind = :string),
                Flopsy.ParameterSpec(:palioxis_root, false, "",
                    "Optional Palioxis root for init"; kind = :string),
                Flopsy.ParameterSpec(
                    :defect_density, false, 1e-3, "Uniform defect density"; kind = :real),
                Flopsy.ParameterSpec(:defect_profile_csv, false, "",
                    "Optional two-column CSV profile (x, defect density) overriding defect_density";
                    kind = :string),
                Flopsy.ParameterSpec(:defect_profile_scale, false, 1.0,
                    "Multiplier applied to values loaded from defect_profile_csv"; kind = :real)
            ],
            "Palioxis-backed trapping backend registered by the Palioxis extension.",
            (data, ctx, reg, block) -> nothing,
            (data,
                ctx,
                reg,
                block) -> begin
                root = String(get(data, "palioxis_root", ""))
                !isempty(root) && Palioxis.init(root)
                pal_model = Palioxis.MultipleDefectModel(String(data["xml_file"]))
                n_gas = pal_model.n_gas
                n_ne = pal_model.n_ne_species
                species = SpeciesInfo[]
                for (i, name) in enumerate(Symbol.(Palioxis.gas_names(pal_model)))
                    push!(species,
                        Flopsy.SpeciesInfo(
                            name, :mobile, :diffusive, i, "Palioxis gas species", true))
                end
                for (j, name) in enumerate(Symbol.(Palioxis.trap_names(pal_model)))
                    push!(species,
                        Flopsy.SpeciesInfo(name, :trapped, :stationary, n_gas + j,
                            "Palioxis trapped species", false))
                end

                build_model = function (mesh, bcs, temperature = nothing)
                    nx = length(mesh.x)
                    defects = _build_defect_profile(data, mesh, pal_model.n_traps)
                    trap_occupancies = _max_trap_occupancies(pal_model, n_ne)
                    trap_groups = Vector{Vector{Int}}()
                    pos = 1
                    for occ in trap_occupancies
                        push!(trap_groups, collect(pos:(pos + occ - 1)))
                        pos += occ
                    end
                    adaptor = Flopsy.HotgatesTrappingAdaptor(
                        collect(1:n_gas),
                        collect((n_gas + 1):(n_gas + n_ne)),
                        Palioxis.gas_names(pal_model),
                        Palioxis.trap_names(pal_model),
                        Palioxis.defect_names(pal_model),
                        Matrix{Float64}(defects),
                        trap_groups
                    )
                    temp = temperature === nothing ? Flopsy.ConstantTemperature(300.0) : temperature
                    diff_coeffs = PalioxisDiffusionCoefficients(pal_model)
                    layout = Flopsy.build_hotgates_variable_layout(adaptor)
                    boundary = Flopsy._build_boundary_operator_from_defs(bcs, layout, diff_coeffs, temp)
                    return Flopsy.build_rd_model(
                        layout = layout,
                        mesh = mesh,
                        reaction = Flopsy.HotgatesReactionOperator(
                            pal_model, adaptor, temp),
                        diffusion = Flopsy.LinearDiffusionOperator(diff_coeffs,
                            layout -> Flopsy.variables_with_tag(layout, :diffusion),
                            nothing, temp),
                        boundary = boundary
                    )
                end

                return Flopsy.BackendDefinition(
                    block.name,
                    :palioxis_trapping,
                    species,
                    build_model,
                    (; palioxis_model = pal_model)
                )
            end,
            :palioxis
        ))

    Flopsy.register_syntax!(registry,
        Flopsy.SyntaxSpec(
            :backend,
            :palioxis_effective_diffusion,
            [
                Flopsy.ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
                Flopsy.ParameterSpec(
                    :xml_file, true, nothing, "Palioxis XML model file"; kind = :string),
                Flopsy.ParameterSpec(:palioxis_root, false, "",
                    "Optional Palioxis root for init"; kind = :string),
                Flopsy.ParameterSpec(
                    :defect_density, false, 1e-3, "Uniform defect density"; kind = :real),
                Flopsy.ParameterSpec(:defect_profile_csv, false, "",
                    "Optional two-column CSV profile (x, defect density) overriding defect_density";
                    kind = :string),
                Flopsy.ParameterSpec(:defect_profile_scale, false, 1.0,
                    "Multiplier applied to values loaded from defect_profile_csv"; kind = :real)
            ],
            "Palioxis-backed equilibrium-trapping effective-diffusion backend.",
            (data, ctx, reg, block) -> nothing,
            (data, ctx, reg, block) -> begin
                root = String(get(data, "palioxis_root", ""))
                !isempty(root) && Palioxis.init(root)
                pal_model = Palioxis.MultipleDefectModel(String(data["xml_file"]))
                n_gas = pal_model.n_gas
                n_gas == 1 || throw(ArgumentError(
                    "palioxis_effective_diffusion currently supports one mobile gas; got $n_gas"))
                species = [
                    Flopsy.SpeciesInfo(Symbol(name), :mobile, :diffusive, i,
                        "Palioxis mobile gas species", true)
                    for (i, name) in enumerate(Palioxis.gas_names(pal_model))
                ]
                build_model = function (mesh, bcs, temperature = nothing)
                    defects = _build_defect_profile(data, mesh, pal_model.n_traps)
                    temp = temperature === nothing ? Flopsy.ConstantTemperature(300.0) : temperature
                    left_bc, right_bc = _mobile_boundary_functions(bcs)
                    return Flopsy.build_palioxis_effective_diffusion_model(
                        palioxis_model = pal_model,
                        mesh = mesh,
                        defects = defects,
                        temperature = temp,
                        left_bc = left_bc,
                        right_bc = right_bc
                    )
                end
                return Flopsy.BackendDefinition(
                    block.name,
                    :palioxis_effective_diffusion,
                    species,
                    build_model,
                    (; palioxis_model = pal_model)
                )
            end,
            :palioxis
        ))

    Flopsy.register_syntax!(registry,
        Flopsy.SyntaxSpec(
            :ic,
            :palioxis_equilibrium,
            [
                Flopsy.ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
                Flopsy.ParameterSpec(
                    :backend, true, nothing, "Referenced backend"; kind = :string),
                Flopsy.ParameterSpec(
                    :driving_quantity, true, nothing, "Supports H_total or mobile";
                    kind = :string, allowed_values = ["H_total", "mobile"]),
                Flopsy.ParameterSpec(:value, true, nothing, "Driving value"; kind = :real),
                Flopsy.ParameterSpec(
                    :temperature, true, nothing, "Equilibrium temperature"; kind = :real)
            ],
            "Palioxis equilibrium initial condition registered by the Palioxis extension.",
            (data,
                ctx,
                reg,
                block) -> begin
                Symbol(data["backend"]) in keys(ctx.backends) ||
                    throw(Flopsy.ConfigValidationError(
                        "Block [ic.$(block.name)] field `backend` references unknown backend `$(data["backend"])`",
                    ))
            end,
            (data,
                ctx,
                reg,
                block) -> begin
                backend = ctx.backends[Symbol(data["backend"])]
                pal_model = backend.metadata.palioxis_model
                affects = Symbol[info.name for info in backend.species]
                return Flopsy.ICDefinition(
                    block.name,
                    block.type_name,
                    Symbol(data["backend"]),
                    affects,
                    (u0,
                        backend_def,
                        model) -> begin
                        value = fill(Float64(data["value"]), model.context.nx)
                        tmp = if data["driving_quantity"] == "H_total"
                            Flopsy.build_ic_from_total_hydrogen(
                                pal_model,
                                model,
                                value,
                                Float64(data["temperature"])
                            )
                        else
                            Flopsy.build_equilibrium_ic(
                                pal_model,
                                model,
                                value,
                                Float64(data["temperature"])
                            )
                        end
                        copyto!(u0, tmp)
                        return u0
                    end
                )
            end,
            :palioxis
        ))

    return registry
end

function _build_defect_profile(data, mesh::Flopsy.Mesh1D, n_traps::Integer)
    profile_csv = String(get(data, "defect_profile_csv", ""))
    if isempty(profile_csv)
        return fill(Float64(data["defect_density"]), n_traps, length(mesh.x))
    end

    table = CSV.File(profile_csv; header = false) |> DataFrames.DataFrame
    DataFrames.ncol(table) >= 2 ||
        throw(ArgumentError("defect_profile_csv must contain at least two columns: x, defect density"))

    xp = Float64.(table[!, 1])
    yp = Float64.(table[!, 2]) .* Float64(get(data, "defect_profile_scale", 1.0))
    length(xp) >= 2 ||
        throw(ArgumentError("defect_profile_csv must contain at least two rows"))
    all(diff(xp) .>= 0) ||
        throw(ArgumentError("defect_profile_csv x values must be sorted in ascending order"))

    vals = [_interp_profile(xp, yp, x) for x in mesh.x]
    return repeat(reshape(vals, 1, :), n_traps, 1)
end

function _interp_profile(xp::AbstractVector, yp::AbstractVector, x::Real)
    x <= xp[1] && return yp[1]
    x >= xp[end] && return yp[end]
    i = searchsortedlast(xp, x)
    x0 = xp[i]
    x1 = xp[i + 1]
    y0 = yp[i]
    y1 = yp[i + 1]
    return y0 + (y1 - y0) * (Float64(x) - x0) / (x1 - x0)
end

function _mobile_boundary_functions(bcs)
    left = nothing
    right = nothing
    for bc in bcs
        bc.method == :weak || throw(ArgumentError(
            "palioxis_effective_diffusion currently supports weak Dirichlet BCs only"))
        fn = let value = bc.value
            t -> value
        end
        if bc.boundary == :left
            left = fn
        elseif bc.boundary == :right
            right = fn
        end
    end
    return left, right
end

function __init__()
    Flopsy.register_plugin_provider!(:palioxis, register_flopsy_plugin!)
end

# ---------------------------------------------------------------------------
# PalioxisDiffusionCoefficients
# ---------------------------------------------------------------------------

"""
    PalioxisDiffusionCoefficients(model)

`AbstractDiffusionCoefficients` backed by a live `Palioxis.MultipleDefectModel`.

Calls `Palioxis.diffusion_constants(model, T)` at every evaluation, so D is
fully consistent with the Palioxis model parameters and the current temperature.
No pre-evaluation at a reference temperature; no approximation.
"""
struct PalioxisDiffusionCoefficients <: Flopsy.AbstractDiffusionCoefficients
    model::Palioxis.MultipleDefectModel
end

function Flopsy.get_D(c::PalioxisDiffusionCoefficients, ivar::Int, ix::Int, T::Real)
    return Palioxis.diffusion_constants(c.model, T)[ivar]
end

"""
    retention_into(out, model, mobile, defects, trapped, T; location=0, thickness=1.0)

Thin Palioxis wrapper used by the equilibrium constitutive evaluator.  It keeps
all retention calls behind a named Flopsy extension-layer function even though
the Julia Palioxis binding exposes the in-place operation as `retention!`.
"""
function retention_into(out::AbstractVector{Float64},
        model::Palioxis.MultipleDefectModel,
        mobile::AbstractVector,
        defects::AbstractVector,
        trapped::AbstractVector,
        T::Real;
        location::Integer = 0,
        thickness::Real = 1.0)
    return Palioxis.retention!(out, model, mobile, defects, trapped, T, location, thickness)
end

"""
    trapped_retention_by_occupation_into(out, model, trap_index, mobile, defect_value, trapped, T)

Thin Palioxis wrapper used by the equilibrium constitutive evaluator.  The
output is ordered as gas-major by occupation for the requested trap, matching
the Palioxis binding contract.
"""
function trapped_retention_by_occupation_into(out::AbstractVector{Float64},
        model::Palioxis.MultipleDefectModel,
        trap_index::Integer,
        mobile::AbstractVector,
        defect_value::Real,
        trapped::AbstractVector,
        T::Real)
    return Palioxis.trapped_retention_by_occupation!(
        out, model, trap_index, mobile, defect_value, trapped, T)
end

"""
    evaluate_equilibrium_state(pal_model, mobile, defects, T) -> PalioxisEquilibriumState

Evaluate the node-local Palioxis equilibrium trapping state for a mobile
concentration vector, local defect vector, and temperature.

The canonical trapped state is obtained from `Palioxis.set_initial_conditions`.
The evaluator then calls the extension-layer `retention_into` and
`trapped_retention_by_occupation_into` wrappers.  `Deff` is computed as
`D_mobile / (d retention_total / d mobile)` using a forward finite difference
of the total-retention constitutive law.  This keeps the solve state mobile-only
while making the diffusion law state dependent.  No state-dependent mass matrix
is introduced in this implementation.
"""
function Flopsy.evaluate_equilibrium_state(
        pal_model::Palioxis.MultipleDefectModel,
        mobile::AbstractVector,
        defects::AbstractVector,
        T::Real)
    mobile_vec = collect(Float64, mobile)
    defects_vec = collect(Float64, defects)
    trapped = Palioxis.set_initial_conditions(pal_model, mobile_vec, T)

    retention_total = zeros(Float64, pal_model.n_gas)
    retention_into(retention_total, pal_model, mobile_vec, defects_vec, trapped, T)

    retention_by_occ = Float64[]
    for trap_index in 0:(pal_model.n_traps - 1)
        max_occ = Int(Palioxis.get_max_trap_occupancy(pal_model, trap_index))
        out = zeros(Float64, pal_model.n_gas * max_occ)
        trapped_retention_by_occupation_into(
            out, pal_model, trap_index, mobile_vec, defects_vec[trap_index + 1],
            trapped, T)
        append!(retention_by_occ, out)
    end

    D_mobile = Palioxis.diffusion_constants(pal_model, T)[1]
    eps_c = max(abs(mobile_vec[1]) * sqrt(eps(Float64)), 1e-30)
    mobile_perturbed = copy(mobile_vec)
    mobile_perturbed[1] += eps_c
    trapped_perturbed = Palioxis.set_initial_conditions(pal_model, mobile_perturbed, T)
    retention_perturbed = zeros(Float64, pal_model.n_gas)
    retention_into(
        retention_perturbed, pal_model, mobile_perturbed, defects_vec,
        trapped_perturbed, T)
    dtotal_dc = (retention_perturbed[1] - retention_total[1]) / eps_c
    Deff = D_mobile / max(dtotal_dc, 1.0)

    return Flopsy.PalioxisEquilibriumState(
        trapped,
        retention_total,
        retention_by_occ,
        Deff
    )
end

# ---------------------------------------------------------------------------
# hotgates_rates! — dispatch for the real Palioxis backend
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# supports_jacobian — enabled for the Palioxis backend
# ---------------------------------------------------------------------------

function Flopsy.supports_jacobian(op::Flopsy.HotgatesReactionOperator{<:Palioxis.MultipleDefectModel})
    return true
end

# ---------------------------------------------------------------------------
# hotgates_jacobian! — analytic per-node Jacobian via Palioxis C API
# ---------------------------------------------------------------------------

function Flopsy.hotgates_jacobian!(
        J_local::AbstractMatrix,
        model::Palioxis.MultipleDefectModel,
        adaptor::Flopsy.HotgatesTrappingAdaptor,
        mobile::AbstractVector,
        defects::AbstractVector,
        trapped::AbstractVector,
        T::Real
)
    jac = Palioxis.time_derivatives_jacobian(model, defects, mobile, trapped, T)

    n_gas = model.n_gas
    n_ne = model.n_ne_species

    # Unpack blocks — stored Fortran col-major, which matches Julia col-major,
    # so reshape gives the correct matrix directly.
    ddx_dx = reshape(jac.ddx_dx.data, n_gas, n_gas)    # d(mobile_rates)/d(mobile)
    ddx_dyc = reshape(jac.ddx_dyc.data, n_gas, n_ne)     # d(mobile_rates)/d(trapped)
    ddydt_dx = reshape(jac.ddydt_dx.data, n_ne, n_gas)    # d(trapped_rates)/d(mobile)
    ddydt_dyc = reshape(jac.ddydt_dyc.data, n_ne, n_ne)     # d(trapped_rates)/d(trapped)

    m_idxs = adaptor.mobile_indices   # row/col positions in local state
    t_idxs = adaptor.trap_indices

    # Accumulate into J_local (row = output variable, col = input variable)
    @inbounds for (im2, cm) in enumerate(m_idxs), (im1, rm) in enumerate(m_idxs)

        J_local[rm, cm] += ddx_dx[im1, im2]
    end
    @inbounds for (it2, ct) in enumerate(t_idxs), (im1, rm) in enumerate(m_idxs)

        J_local[rm, ct] += ddx_dyc[im1, it2]
    end
    @inbounds for (im2, cm) in enumerate(m_idxs), (it1, rt) in enumerate(t_idxs)

        J_local[rt, cm] += ddydt_dx[it1, im2]
    end
    @inbounds for (it2, ct) in enumerate(t_idxs), (it1, rt) in enumerate(t_idxs)

        J_local[rt, ct] += ddydt_dyc[it1, it2]
    end

    return J_local
end

# ---------------------------------------------------------------------------
# hotgates_rates!
# ---------------------------------------------------------------------------

function Flopsy.hotgates_rates!(
        dmobile::AbstractVector,
        dtrapped::AbstractVector,
        model::Palioxis.MultipleDefectModel,
        mobile::AbstractVector,
        defects::AbstractVector,
        trapped::AbstractVector,
        T::Real
)
    Palioxis.rates!(dmobile, dtrapped, model, mobile, defects, trapped, T)
    return nothing
end

# ---------------------------------------------------------------------------
# build_palioxis_effective_diffusion_model
# ---------------------------------------------------------------------------

"""
    build_palioxis_effective_diffusion_model(; palioxis_model, mesh, defects, temperature,
                                               left_bc=nothing, right_bc=nothing)

Build a mobile-only effective-diffusion model backed by Palioxis equilibrium
trapping.  The solve state contains only Palioxis gas/mobile variables tagged
with `:diffusion`; trapped populations are evaluated as auxiliary equilibrium
quantities through `evaluate_equilibrium_state`.

This builder deliberately does not mix dynamic trapped variables into the
effective-diffusion model family.  Use `build_dynamic_ic_from_stage1` to hand a
completed effective-diffusion result into `build_palioxis_trapping_model`.
"""
function Flopsy.build_palioxis_effective_diffusion_model(;
        palioxis_model::Palioxis.MultipleDefectModel,
        mesh::Flopsy.Mesh1D,
        defects::AbstractMatrix,
        temperature::Flopsy.AbstractTemperatureProvider,
        left_bc = nothing,
        right_bc = nothing)
    palioxis_model.n_gas == 1 || throw(ArgumentError(
        "build_palioxis_effective_diffusion_model currently supports one mobile gas"))
    nx = length(mesh.x)
    size(defects) == (palioxis_model.n_traps, nx) ||
        throw(DimensionMismatch(
            "defects must have shape ($(palioxis_model.n_traps), $nx), got $(size(defects))"))

    vars = [
        Flopsy.VariableInfo(Symbol(name), :mobile, Set([:diffusion]))
        for name in Palioxis.gas_names(palioxis_model)
    ]
    layout = Flopsy.VariableLayout(vars)
    selector = layout -> Flopsy.variables_with_tag(layout, :diffusion)
    evaluator = (mobile, defects, T) -> Flopsy.evaluate_equilibrium_state(
        palioxis_model, mobile, defects, T)
    diffusion = Flopsy.NonlinearDiffusionOperator(
        selector, evaluator, Matrix{Float64}(defects), temperature;
        left = left_bc, right = right_bc)

    return Flopsy.build_rd_model(
        layout = layout,
        mesh = mesh,
        diffusion = diffusion,
        aux = (
            model_family = :palioxis_effective_diffusion,
            palioxis_model = palioxis_model,
            defects = Matrix{Float64}(defects),
            temperature = temperature,
            gas_names = Palioxis.gas_names(palioxis_model),
            trap_names = Palioxis.trap_names(palioxis_model),
            defect_names = Palioxis.defect_names(palioxis_model),
            occupancy_ordering = _occupancy_ordering(palioxis_model)
        )
    )
end

# ---------------------------------------------------------------------------
# build_palioxis_trapping_model
# ---------------------------------------------------------------------------

function Flopsy.build_palioxis_trapping_model(;
        palioxis_model::Palioxis.MultipleDefectModel,
        mesh::Flopsy.Mesh1D,
        defects::AbstractMatrix,
        temperature::Flopsy.AbstractTemperatureProvider,
        left_bc = nothing,
        right_bc = nothing
)
    nx = length(mesh.x)

    size(defects) == (palioxis_model.n_traps, nx) ||
        throw(DimensionMismatch(
            "defects must have shape ($(palioxis_model.n_traps), $nx), " *
            "got $(size(defects))"
        ))

    n_gas = palioxis_model.n_gas
    n_ne = palioxis_model.n_ne_species

    # Build trap_groups from the Palioxis occupancy structure.
    # Each defect type has max_occ[d] fill levels; levels within a type are coupled
    # (multi-occupancy trapping couples adjacent occupancy states).
    trap_occupancies = _max_trap_occupancies(palioxis_model, n_ne)
    trap_groups = Vector{Vector{Int}}()
    pos = 1
    for occ in trap_occupancies
        push!(trap_groups, collect(pos:(pos + occ - 1)))
        pos += occ
    end

    adaptor = Flopsy.HotgatesTrappingAdaptor(
        collect(1:n_gas),
        collect((n_gas + 1):(n_gas + n_ne)),
        Palioxis.gas_names(palioxis_model),
        Palioxis.trap_names(palioxis_model),
        Palioxis.defect_names(palioxis_model),
        Matrix{Float64}(defects),
        trap_groups
    )

    # Diffusion coefficients come from Palioxis at each time step — fully T-dependent.
    # The coefficients object wraps the Palioxis model; the selector restricts
    # application to mobile (diffusing) variables only.
    diff_coeffs = PalioxisDiffusionCoefficients(palioxis_model)

    # Flat coefficient vector for LinearDiffusionOperator: the selector only
    # accesses mobile indices so trapped slots are never read, but we need
    # the *diffusion* operator to store the same PalioxisDiffusionCoefficients.
    selector = layout -> Flopsy.variables_with_tag(layout, :diffusion)

    diffusion = Flopsy.LinearDiffusionOperator(diff_coeffs, selector, nothing, temperature)

    boundary = if left_bc !== nothing || right_bc !== nothing
        Flopsy.WeakDirichletBoundaryOperator(selector, diff_coeffs, temperature;
            left = left_bc, right = right_bc)
    else
        nothing
    end

    return Flopsy.build_rd_model(
        layout = Flopsy.build_hotgates_variable_layout(adaptor),
        mesh = mesh,
        reaction = Flopsy.HotgatesReactionOperator(palioxis_model, adaptor, temperature),
        diffusion = diffusion,
        boundary = boundary
    )
end

function _max_trap_occupancies(palioxis_model::Palioxis.MultipleDefectModel, n_ne::Integer)
    n_ne == 0 && return Int[]
    return [Palioxis.get_max_trap_occupancy(palioxis_model, i)
            for i in 0:(palioxis_model.n_traps - 1)]
end

function _occupancy_ordering(palioxis_model::Palioxis.MultipleDefectModel)
    labels = String[]
    for trap_index in 0:(palioxis_model.n_traps - 1)
        max_occ = Int(Palioxis.get_max_trap_occupancy(palioxis_model, trap_index))
        for gas in Palioxis.gas_names(palioxis_model), occ in 1:max_occ
            push!(labels, "trap=$(trap_index);gas=$(gas);occupation=$(occ)")
        end
    end
    return labels
end

function _effective_aux(result::Flopsy.SimulationResult)
    aux = result.model.context.aux
    if aux isa NamedTuple && hasproperty(aux, :model_family) &&
       getproperty(aux, :model_family) == :palioxis_effective_diffusion
        return aux
    end
    throw(ArgumentError("result is not a Palioxis effective-diffusion result"))
end

"""
    compute_equilibrium_aux_fields(result, pal_model=nothing)

Compute pointwise auxiliary equilibrium fields for a Palioxis effective-diffusion
result.  Arrays are ordered `(time, node, component)` except `Deff`, which is
`(time, node)`.  Metadata includes a schema version, trap indices, isotope
ordering, occupation ordering, and the saved temperature field.
"""
function Flopsy.compute_equilibrium_aux_fields(
        result::Flopsy.SimulationResult,
        pal_model::Palioxis.MultipleDefectModel)
    aux = _effective_aux(result)
    pal = pal_model
    defects = aux.defects
    temperature = aux.temperature
    layout = result.model.layout
    nx = result.model.context.nx
    nt = length(result.solution.u)
    n_gas = pal.n_gas
    n_ne = pal.n_ne_species
    n_occ = length(aux.occupancy_ordering)

    trapped = zeros(Float64, nt, nx, n_ne)
    retention_total = zeros(Float64, nt, nx, n_gas)
    retention_by_occ = zeros(Float64, nt, nx, n_occ)
    Deff = zeros(Float64, nt, nx)
    temperature_field = zeros(Float64, nt, nx)
    vars = Flopsy.variables_with_tag(layout, :diffusion)

    for it in 1:nt
        t = result.solution.t[it]
        U = Flopsy.state_view(result.solution.u[it], layout, nx)
        for ix in 1:nx
            T = Float64(Flopsy.temperature_at(temperature, result.model.context, t, ix))
            mobile = [Float64(U[ivar, ix]) for ivar in vars]
            eq = Flopsy.evaluate_equilibrium_state(pal, mobile, @view(defects[:, ix]), T)
            trapped[it, ix, :] .= eq.trapped
            retention_total[it, ix, :] .= eq.retention_total
            retention_by_occ[it, ix, :] .= eq.retention_by_occupation
            Deff[it, ix] = eq.Deff
            temperature_field[it, ix] = T
        end
    end

    metadata = Dict{String, Any}(
        "schema_version" => "flopsy.equilibrium_aux.v1",
        "trap_indices" => join(0:(pal.n_traps - 1), ","),
        "isotope_ordering" => join(Palioxis.gas_names(pal), ","),
        "trap_ordering" => join(Palioxis.trap_names(pal), ","),
        "occupancy_ordering" => join(aux.occupancy_ordering, "|"),
        "temperature_field" => "temperature_K"
    )

    return Dict{Symbol, Any}(
        :equilibrium_trapped => trapped,
        :retention_total => retention_total,
        :retention_by_occupation => retention_by_occ,
        :Deff => Deff,
        :temperature_K => temperature_field,
        :metadata => metadata
    )
end

# ---------------------------------------------------------------------------
# build_equilibrium_ic — mobile profile → Palioxis equilibrium trapped
# ---------------------------------------------------------------------------

function Flopsy.build_equilibrium_ic(
        palioxis_model::Palioxis.MultipleDefectModel,
        model::Flopsy.SystemModel,
        mobile_profile::AbstractVecOrMat,
        T::Real
)
    adaptor = model.operators.reaction.adaptor
    n_gas = palioxis_model.n_gas
    n_ne = palioxis_model.n_ne_species
    nvars = n_gas + n_ne
    nx = size(adaptor.defects, 2)

    mobile_mat = if mobile_profile isa AbstractVector
        n_gas == 1 || throw(ArgumentError(
            "Vector mobile_profile only valid for single-gas models (n_gas=$(n_gas))"
        ))
        reshape(Float64.(mobile_profile), 1, nx)
    else
        Matrix{Float64}(mobile_profile)
    end

    size(mobile_mat) == (n_gas, nx) || throw(DimensionMismatch(
        "mobile_profile must have shape ($n_gas, $nx), got $(size(mobile_mat))"
    ))

    u0 = zeros(Float64, nvars * nx)
    U0 = reshape(u0, nvars, nx)

    for ix in 1:nx
        mobile = mobile_mat[:, ix]
        trapped = Palioxis.set_initial_conditions(palioxis_model, mobile, T)

        for (j, idx) in enumerate(adaptor.mobile_indices)
            U0[idx, ix] = mobile[j]
        end
        for (j, idx) in enumerate(adaptor.trap_indices)
            U0[idx, ix] = trapped[j]
        end
    end

    return u0
end

# ---------------------------------------------------------------------------
# build_ic_from_total_hydrogen — total H profile → Palioxis steady-state
# ---------------------------------------------------------------------------

function Flopsy.build_ic_from_total_hydrogen(
        palioxis_model::Palioxis.MultipleDefectModel,
        model::Flopsy.SystemModel,
        total_hydrogen::AbstractVector,
        T::Real
)
    adaptor = model.operators.reaction.adaptor
    n_gas = palioxis_model.n_gas
    n_ne = palioxis_model.n_ne_species
    nvars = n_gas + n_ne
    nx = size(adaptor.defects, 2)

    length(total_hydrogen) == nx || throw(DimensionMismatch(
        "total_hydrogen must have length nx=$nx, got $(length(total_hydrogen))"
    ))

    u0 = zeros(Float64, nvars * nx)
    U0 = reshape(u0, nvars, nx)

    for ix in 1:nx
        C_tot = Float64(total_hydrogen[ix])
        defects = adaptor.defects[:, ix]

        # Initial guess: divide evenly across all DOFs, then clamp to valid bounds.
        frac = C_tot / (n_gas + n_ne)
        mobile = fill(frac, n_gas)
        trapped = fill(frac, n_ne)

        Palioxis.ensure_bounds!(palioxis_model, mobile, collect(Float64, defects), trapped)

        # Iterate to the nearest equilibrium point (conserves chemistry, not total H).
        Palioxis.calculate_steady_state!(palioxis_model, mobile, defects, trapped, T)

        for (j, idx) in enumerate(adaptor.mobile_indices)
            U0[idx, ix] = mobile[j]
        end
        for (j, idx) in enumerate(adaptor.trap_indices)
            U0[idx, ix] = trapped[j]
        end
    end

    return u0
end

function _final_temperature(result::Flopsy.SimulationResult)
    aux = _effective_aux(result)
    t = result.solution.t[end]
    return Float64(Flopsy.temperature_at(aux.temperature, result.model.context, t, 1))
end

"""
    build_dynamic_ic_from_stage1(result; mode=:recompute_from_mobile)

Build the initial-condition vector for the existing dynamic Palioxis model from
a Palioxis effective-diffusion result.  The returned vector uses the dynamic
state ordering expected by `build_palioxis_trapping_model`: all mobile gas
variables followed by all dynamic trapped variables at each node.

Supported modes:
- `:use_exported_trapped` — use `result.summaries[:equilibrium_aux]` trapped
  fields directly; no Palioxis recomputation is performed.
- `:recompute_from_mobile` — recompute equilibrium trapped values from the final
  mobile profile and final temperature.
- `:use_exported_retention` — reserved for future reconstruction rules; throws
  unless such a rule is supplied by a later implementation.
"""
function Flopsy.build_dynamic_ic_from_stage1(
        result::Flopsy.SimulationResult; mode::Symbol = :recompute_from_mobile)
    aux = _effective_aux(result)
    pal = aux.palioxis_model
    nx = result.model.context.nx
    n_gas = pal.n_gas
    n_ne = pal.n_ne_species
    nvars = n_gas + n_ne
    U_mobile = Flopsy.state_view(result.solution.u[end], result.model.layout, nx)
    mobile_vars = Flopsy.variables_with_tag(result.model.layout, :diffusion)

    u0 = zeros(Float64, nvars * nx)
    U0 = reshape(u0, nvars, nx)

    trapped_source = nothing
    if mode == :use_exported_trapped
        haskey(result.summaries, :equilibrium_aux) ||
            throw(ArgumentError(":use_exported_trapped requires result.summaries[:equilibrium_aux]"))
        eq_aux = result.summaries[:equilibrium_aux]
        trapped_source = eq_aux[:equilibrium_trapped]
    elseif mode == :use_exported_retention
        throw(ArgumentError(
            ":use_exported_retention is only valid when a retention reconstruction rule is defined"))
    elseif mode != :recompute_from_mobile
        throw(ArgumentError("Unknown stage handoff mode :$mode"))
    end

    T = _final_temperature(result)
    for ix in 1:nx
        mobile = [Float64(U_mobile[ivar, ix]) for ivar in mobile_vars]
        U0[1:n_gas, ix] .= mobile
        if mode == :use_exported_trapped
            U0[(n_gas + 1):nvars, ix] .= @view trapped_source[end, ix, :]
        else
            U0[(n_gas + 1):nvars, ix] .= Palioxis.set_initial_conditions(pal, mobile, T)
        end
    end

    return u0
end

"""
    run_effective_diffusion(; palioxis_model, mesh, defects, temperature, u0,
                              tspan, solver_config, left_bc=nothing, right_bc=nothing,
                              compute_aux=true)

Build and run a Palioxis effective-diffusion stage.  When `compute_aux=true`,
equilibrium auxiliary fields are attached to the returned result under
`result.summaries[:equilibrium_aux]`.
"""
function Flopsy.run_effective_diffusion(;
        palioxis_model::Palioxis.MultipleDefectModel,
        mesh::Flopsy.Mesh1D,
        defects::AbstractMatrix,
        temperature::Flopsy.AbstractTemperatureProvider,
        u0,
        tspan,
        solver_config::Flopsy.SolverConfig,
        left_bc = nothing,
        right_bc = nothing,
        compute_aux::Bool = true)
    model = Flopsy.build_palioxis_effective_diffusion_model(
        palioxis_model = palioxis_model,
        mesh = mesh,
        defects = defects,
        temperature = temperature,
        left_bc = left_bc,
        right_bc = right_bc)
    sol = Flopsy.solve_problem(model, u0, tspan, solver_config)
    result = Flopsy.wrap_result(model, sol, solver_config;
        metadata = Dict{String, Any}("model_type" => "palioxis_effective_diffusion"))
    compute_aux && Flopsy.attach_equilibrium_aux_fields!(result)
    return result
end

"""
    run_dynamic_palioxis(; palioxis_model, mesh, defects, temperature, u0, tspan,
                          solver_config, left_bc=nothing, right_bc=nothing)

Build and run the existing dynamic Palioxis trapping model from an explicit
initial-condition vector.
"""
function Flopsy.run_dynamic_palioxis(;
        palioxis_model::Palioxis.MultipleDefectModel,
        mesh::Flopsy.Mesh1D,
        defects::AbstractMatrix,
        temperature::Flopsy.AbstractTemperatureProvider,
        u0,
        tspan,
        solver_config::Flopsy.SolverConfig,
        left_bc = nothing,
        right_bc = nothing)
    model = Flopsy.build_palioxis_trapping_model(
        palioxis_model = palioxis_model,
        mesh = mesh,
        defects = defects,
        temperature = temperature,
        left_bc = left_bc,
        right_bc = right_bc)
    sol = Flopsy.solve_problem(model, u0, tspan, solver_config)
    return Flopsy.wrap_result(model, sol, solver_config;
        metadata = Dict{String, Any}("model_type" => "palioxis_dynamic"))
end

"""
    run_implantation_then_desorption(; kwargs...)

Run the two-stage Palioxis workflow:
1. mobile-only effective diffusion with equilibrium trapping auxiliaries,
2. dynamic Palioxis trapping/desorption initialised from stage one.

Returns `(implantation=result1, desorption=result2, dynamic_ic=u0_dynamic)`.
"""
function Flopsy.run_implantation_then_desorption(;
        palioxis_model::Palioxis.MultipleDefectModel,
        mesh::Flopsy.Mesh1D,
        defects::AbstractMatrix,
        implantation_temperature::Flopsy.AbstractTemperatureProvider,
        desorption_temperature::Flopsy.AbstractTemperatureProvider,
        implantation_u0,
        implantation_tspan,
        desorption_tspan,
        implantation_solver_config::Flopsy.SolverConfig,
        desorption_solver_config::Flopsy.SolverConfig,
        implantation_left_bc = nothing,
        implantation_right_bc = nothing,
        desorption_left_bc = nothing,
        desorption_right_bc = nothing,
        handoff_mode::Symbol = :use_exported_trapped)
    result1 = Flopsy.run_effective_diffusion(
        palioxis_model = palioxis_model,
        mesh = mesh,
        defects = defects,
        temperature = implantation_temperature,
        u0 = implantation_u0,
        tspan = implantation_tspan,
        solver_config = implantation_solver_config,
        left_bc = implantation_left_bc,
        right_bc = implantation_right_bc,
        compute_aux = true)
    u0_dynamic = Flopsy.build_dynamic_ic_from_stage1(result1; mode = handoff_mode)
    result2 = Flopsy.run_dynamic_palioxis(
        palioxis_model = palioxis_model,
        mesh = mesh,
        defects = defects,
        temperature = desorption_temperature,
        u0 = u0_dynamic,
        tspan = desorption_tspan,
        solver_config = desorption_solver_config,
        left_bc = desorption_left_bc,
        right_bc = desorption_right_bc)
    return (implantation = result1, desorption = result2, dynamic_ic = u0_dynamic)
end

end # module PalioxisExt
