"""
    PalioxisEquilibriumState

Node-local equilibrium-trapping constitutive state used by the Palioxis-backed
effective-diffusion model.

Fields:
- `trapped` — equilibrium non-equilibrium-trap state vector returned by Palioxis
- `retention_total` — total retention by mobile isotope/gas
- `retention_by_occupation` — occupation-resolved trapped retention, flattened in
  trap-major, gas-major, occupation-major order
- `Deff` — scalar effective diffusivity used by `NonlinearDiffusionOperator`
"""
struct PalioxisEquilibriumState
    trapped::Vector{Float64}
    retention_total::Vector{Float64}
    retention_by_occupation::Vector{Float64}
    Deff::Float64
end

"""
    build_palioxis_effective_diffusion_model(; kwargs...)

Build a Palioxis-backed effective-diffusion model for equilibrium trapping.

This is a package-extension entry point implemented by `ext/PalioxisExt.jl`.
The resulting model contains only mobile variables in the solve state; trapped
populations and retention quantities are evaluated as auxiliary equilibrium
fields.
"""
function build_palioxis_effective_diffusion_model end

"""
    evaluate_equilibrium_state(pal_model, mobile, defects, T)

Evaluate a backend-specific local equilibrium trapping state.

The Palioxis extension implements this method and returns a
`PalioxisEquilibriumState` containing equilibrium trapped concentrations,
retention summaries, and the scalar effective diffusivity used by the nonlinear
diffusion operator.
"""
function evaluate_equilibrium_state end

"""
    compute_equilibrium_aux_fields(result, pal_model=nothing)

Compute saved-time pointwise auxiliary fields for an equilibrium-trapping result.

The default implementation returns an empty dictionary.  Backend extensions
override this for models that carry enough metadata to reconstruct equilibrium
trapped populations, retention, and effective diffusivity.
"""
function compute_equilibrium_aux_fields(result::SimulationResult, pal_model = nothing)
    if pal_model !== nothing
        return Dict{Symbol, Any}()
    end
    aux = result.model.context.aux
    if aux isa NamedTuple && hasproperty(aux, :palioxis_model)
        return compute_equilibrium_aux_fields(result, getproperty(aux, :palioxis_model))
    end
    return Dict{Symbol, Any}()
end

"""
    attach_equilibrium_aux_fields!(result, pal_model=nothing) -> result

Compute equilibrium auxiliary fields and store them under
`result.summaries[:equilibrium_aux]`.
"""
function attach_equilibrium_aux_fields!(result::SimulationResult, pal_model = nothing)
    result.summaries[:equilibrium_aux] = compute_equilibrium_aux_fields(result, pal_model)
    return result
end

"""
    build_dynamic_ic_from_stage1(result; mode=:recompute_from_mobile)

Construct a dynamic Palioxis initial condition from an effective-diffusion
stage-one result.  Implemented by the Palioxis extension.
"""
function build_dynamic_ic_from_stage1 end

"""
    run_effective_diffusion(; kwargs...)

Convenience helper for running a Palioxis effective-diffusion stage.

Implemented by the Palioxis extension.
"""
function run_effective_diffusion end

"""
    run_dynamic_palioxis(; kwargs...)

Convenience helper for running the dynamic Palioxis model, usually after
`build_dynamic_ic_from_stage1`. Implemented by the Palioxis extension.
"""
function run_dynamic_palioxis end

"""
    run_implantation_then_desorption(; kwargs...)

Run a two-stage workflow: effective-diffusion implantation/loading followed by
dynamic Palioxis desorption.  Implemented by the Palioxis extension.
"""
function run_implantation_then_desorption end
