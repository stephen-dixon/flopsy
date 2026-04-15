function register_builtin_syntax!(registry::SyntaxRegistry)
    _register_mesh_syntax!(registry)
    _register_backend_syntax!(registry)
    _register_ic_syntax!(registry)
    _register_bc_syntax!(registry)
    _register_output_syntax!(registry)
    _register_problem_syntax!(registry)
    return registry
end

function _register_mesh_syntax!(registry::SyntaxRegistry)
    register_syntax!(registry, SyntaxSpec(
        :mesh,
        :uniform_1d,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:xmin, true, nothing, "Domain minimum"; kind = :real),
            ParameterSpec(:xmax, true, nothing, "Domain maximum"; kind = :real),
            ParameterSpec(:nx, true, nothing, "Number of nodes"; kind = :integer),
        ],
        "Uniform one-dimensional mesh.",
        (data, ctx, reg, block) -> begin
            Float64(data["xmax"]) > Float64(data["xmin"]) ||
                throw(ConfigValidationError("Block [mesh.$(block.name)] field `xmax` must be greater than `xmin`"))
            Int(data["nx"]) >= 2 || throw(ConfigValidationError("Block [mesh.$(block.name)] field `nx` must be at least 2"))
        end,
        (data, ctx, reg, block) -> Mesh1D(Float64(data["xmax"]) - Float64(data["xmin"]), Int(data["nx"])),
        :builtin,
    ))

    register_syntax!(registry, SyntaxSpec(
        :mesh,
        :biased_1d,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:xmin, true, nothing, "Domain minimum"; kind = :real),
            ParameterSpec(:xmax, true, nothing, "Domain maximum"; kind = :real),
            ParameterSpec(:nx, true, nothing, "Number of nodes"; kind = :integer),
            ParameterSpec(:bias_kind, true, nothing, "Bias strategy"; kind = :string),
            ParameterSpec(:ratio, true, nothing, "Bias ratio"; kind = :real),
            ParameterSpec(:cluster, true, nothing, "Clustering location"; kind = :string),
        ],
        "Experimental biased 1D mesh syntax. Numerics are not implemented yet.",
        (data, ctx, reg, block) -> throw(ConfigValidationError(
            "Block [mesh.$(block.name)] type `biased_1d` is experimental and not available for production use yet",
        )),
        (data, ctx, reg, block) -> nothing,
        :builtin,
    ))
end

function _register_backend_syntax!(registry::SyntaxRegistry)
    register_syntax!(registry, SyntaxSpec(
        :backend,
        :diffusion_1d,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:diffusion_coefficient, true, nothing, "Diffusion coefficient"; kind = :real),
            ParameterSpec(:reaction_rate, false, 0.0, "Uniform decay rate"; kind = :real),
            ParameterSpec(:species_name, false, "u", "Species name"; kind = :string),
        ],
        "Single-species reaction-diffusion backend with one diffusive field.",
        (data, ctx, reg, block) -> nothing,
        (data, ctx, reg, block) -> _build_diffusion_backend(block.name, data),
        :builtin,
    ))

    register_syntax!(registry, SyntaxSpec(
        :backend,
        :trapping_1d,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:k_trap, true, nothing, "Trap rate"; kind = :real),
            ParameterSpec(:k_detrap, true, nothing, "Detrap rate"; kind = :real),
            ParameterSpec(:diffusion_coefficient, true, nothing, "Mobile diffusion coefficient"; kind = :real),
            ParameterSpec(:mobile_species, false, "H_mobile", "Mobile species name"; kind = :string),
            ParameterSpec(:trapped_species, false, "H_trapped", "Trapped species name"; kind = :string),
        ],
        "Two-species trapping backend with one diffusive mobile species and one stationary trapped species.",
        (data, ctx, reg, block) -> nothing,
        (data, ctx, reg, block) -> _build_trapping_backend(block.name, data),
        :builtin,
    ))

    register_syntax!(registry, SyntaxSpec(
        :backend,
        :hotgates_trapping,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:k_trap, true, nothing, "Trap rate"; kind = :real),
            ParameterSpec(:k_detrap, true, nothing, "Detrap rate"; kind = :real),
            ParameterSpec(:diffusion_coefficient, true, nothing, "Mobile diffusion coefficient"; kind = :real),
            ParameterSpec(:temperature, false, 300.0, "Uniform temperature"; kind = :real),
            ParameterSpec(:mobile_species, false, "H_mobile", "Mobile species name"; kind = :string),
            ParameterSpec(:trapped_species, false, "H_trapped", "Trapped species name"; kind = :string),
        ],
        "Fake Hotgates-style backend for testing and config-driven workflows without Palioxis.",
        (data, ctx, reg, block) -> nothing,
        (data, ctx, reg, block) -> _build_fake_hotgates_backend(block.name, data),
        :builtin,
    ))
end

function _register_ic_syntax!(registry::SyntaxRegistry)
    register_syntax!(registry, SyntaxSpec(
        :ic,
        :uniform_species,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:species, true, nothing, "Target species"; kind = :string),
            ParameterSpec(:value, true, nothing, "Uniform value"; kind = :real),
            ParameterSpec(:backend, false, nothing, "Optional backend name"; kind = :string),
        ],
        "Uniform species-targeting initial condition.",
        (data, ctx, reg, block) -> nothing,
        (data, ctx, reg, block) -> ICDefinition(
            block.name,
            block.type_name,
            haskey(data, "backend") ? Symbol(data["backend"]) : nothing,
            [Symbol(data["species"])],
            (u0, backend, model) -> begin
                ivar = _species_index(backend, Symbol(data["species"]))
                U0 = state_view(u0, model.layout, model.context.nx)
                U0[ivar, :] .= Float64(data["value"])
                return u0
            end,
        ),
        :builtin,
    ))
end

function _register_bc_syntax!(registry::SyntaxRegistry)
    register_syntax!(registry, SyntaxSpec(
        :bc,
        :dirichlet,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:species, true, nothing, "Target species"; kind = :string),
            ParameterSpec(:boundary, true, nothing, "left or right"; kind = :string, allowed_values = ["left", "right"]),
            ParameterSpec(:value, true, nothing, "Boundary value"; kind = :real),
            ParameterSpec(
                :method,
                false,
                "weak",
                "weak, penalty, mass_matrix, callback, eliminated";
                kind = :string,
                allowed_values = ["weak", "penalty", "mass_matrix", "callback", "eliminated"],
            ),
        ],
        "Species-targeted Dirichlet boundary condition. If no BC blocks are provided, the default is implicit zero-flux Neumann behaviour from the diffusion operator.",
        (data, ctx, reg, block) -> nothing,
        (data, ctx, reg, block) -> BCDefinition(
            block.name,
            block.type_name,
            Symbol(data["species"]),
            Symbol(data["boundary"]),
            Float64(data["value"]),
            Symbol(data["method"]),
        ),
        :builtin,
    ))
end

function _register_output_syntax!(registry::SyntaxRegistry)
    register_syntax!(registry, SyntaxSpec(
        :output,
        :hdf5,
        [
            ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
            ParameterSpec(:file, true, nothing, "Output HDF5 file"; kind = :string),
            ParameterSpec(:xdmf, false, false, "Generate XDMF companion"; kind = :bool),
            ParameterSpec(:summary_csv, false, "", "Optional summary CSV path"; kind = :string),
        ],
        "HDF5 field output with optional XDMF companion generation.",
        (data, ctx, reg, block) -> nothing,
        (data, ctx, reg, block) -> OutputDefinition(
            block.name,
            block.type_name,
            String(data["file"]),
            Bool(data["xdmf"]),
        ),
        :builtin,
    ))
end

function _register_problem_syntax!(registry::SyntaxRegistry)
    params = [
        ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
        ParameterSpec(:mesh, true, nothing, "Referenced mesh object"; kind = :string),
        ParameterSpec(:backend, true, nothing, "Referenced backend object"; kind = :string),
        ParameterSpec(:ics, false, String[], "Referenced IC objects"; kind = :vector, element_kind = :string),
        ParameterSpec(:bcs, false, String[], "Referenced BC objects"; kind = :vector, element_kind = :string),
        ParameterSpec(:outputs, false, String[], "Referenced output objects"; kind = :vector, element_kind = :string),
        ParameterSpec(:tspan, true, nothing, "Two-element time span"; kind = :vector, element_kind = :real),
        ParameterSpec(:formulation, false, "unsplit", "Solver formulation"; kind = :string, allowed_values = ["unsplit", "imex", "imex_reaction", "split", "residual"]),
        ParameterSpec(:algorithm, false, "Rodas5", "Solver algorithm"; kind = :string, allowed_values = ["Rodas5", "Rodas5P", "KenCarp4", "CVODE_BDF"]),
        ParameterSpec(:abstol, false, 1e-8, "Absolute tolerance"; kind = :real),
        ParameterSpec(:reltol, false, 1e-6, "Relative tolerance"; kind = :real),
        ParameterSpec(:saveat, false, nothing, "Save times"; kind = :vector, element_kind = :real),
        ParameterSpec(:dt, false, nothing, "Fixed/split time step"; kind = :real),
    ]

    for type_name in (:simulation, :diffusion_1d, :trapping_1d, :hotgates_trapping)
        register_syntax!(registry, SyntaxSpec(
            :problem,
            type_name,
            params,
            "Simulation assembly block referencing named mesh/backend/ic/bc/output objects.",
            (data, ctx, reg, block) -> begin
                length(data["tspan"]) == 2 ||
                    throw(ConfigValidationError("Block [problem.$(block.name)] field `tspan` must contain exactly two values"))
            end,
            (data, ctx, reg, block) -> ProblemDefinition(
                block.name,
                block.type_name,
                Symbol(data["mesh"]),
                Symbol(data["backend"]),
                Symbol.(get(data, "ics", String[])),
                Symbol.(get(data, "bcs", String[])),
                Symbol.(get(data, "outputs", String[])),
                (Float64(data["tspan"][1]), Float64(data["tspan"][2])),
                InputSolverConfig(
                    Symbol(lowercase(String(data["formulation"]))),
                    Symbol(data["algorithm"]),
                    get(data, "dt", nothing) === nothing ? nothing : Float64(data["dt"]),
                    Float64(data["abstol"]),
                    Float64(data["reltol"]),
                    get(data, "saveat", nothing) === nothing ? nothing : Float64.(data["saveat"]),
                ),
            ),
            :builtin,
        ))
    end
end

function _build_diffusion_backend(name::Symbol, data::Dict{String, Any})
    species_name = Symbol(data["species_name"])
    species = [SpeciesInfo(species_name, :field, :diffusive, 1, "Primary diffusive field", true)]

    build_model = function (mesh, bcs)
        layout = VariableLayout([VariableInfo(species_name, :field, Set([:diffusion, :reaction]))])
        selector(layout::VariableLayout) = [1]
        coeffs = ConstantDiffusion([Float64(data["diffusion_coefficient"])])
        reaction = ToyReactionOperator(Float64(data["reaction_rate"]))
        diffusion = LinearDiffusionOperator(coeffs, selector, nothing)
        boundary = _build_boundary_operator_from_defs(bcs, layout, coeffs)
        return build_rd_model(
            layout = layout,
            mesh = mesh,
            reaction = reaction,
            diffusion = diffusion,
            boundary = boundary,
        )
    end

    return BackendDefinition(name, :diffusion_1d, species, build_model, (;))
end

function _build_trapping_backend(name::Symbol, data::Dict{String, Any})
    mobile_name = Symbol(data["mobile_species"])
    trapped_name = Symbol(data["trapped_species"])
    species = [
        SpeciesInfo(mobile_name, :mobile, :diffusive, 1, "Mobile diffusive species", true),
        SpeciesInfo(trapped_name, :trapped, :stationary, 2, "Trapped stationary species", false),
    ]

    build_model = function (mesh, bcs)
        layout = trapping_variable_layout(mobile_name = mobile_name, trap_name = trapped_name)
        reaction = SimpleTrappingReactionOperator(Float64(data["k_trap"]), Float64(data["k_detrap"]), 1, 2)
        coeffs = ConstantDiffusion([Float64(data["diffusion_coefficient"]), 0.0])
        selector(layout::VariableLayout) = variables_with_tag(layout, :diffusion)
        diffusion = LinearDiffusionOperator(coeffs, selector, nothing)
        boundary = _build_boundary_operator_from_defs(bcs, layout, coeffs)
        return build_rd_model(
            layout = layout,
            mesh = mesh,
            reaction = reaction,
            diffusion = diffusion,
            boundary = boundary,
        )
    end

    return BackendDefinition(name, :trapping_1d, species, build_model, (;))
end

function _build_fake_hotgates_backend(name::Symbol, data::Dict{String, Any})
    mobile_name = String(data["mobile_species"])
    trapped_name = String(data["trapped_species"])
    species = [
        SpeciesInfo(Symbol(mobile_name), :mobile, :diffusive, 1, "Mobile diffusive species", true),
        SpeciesInfo(Symbol(trapped_name), :trapped, :stationary, 2, "Trapped stationary species", false),
    ]

    build_model = function (mesh, bcs)
        nx = length(mesh.x)
        adaptor = HotgatesTrappingAdaptor([1], [2], [mobile_name], [trapped_name], String[], zeros(Float64, 0, nx))
        backend = FakeHotgatesModel(Float64(data["k_trap"]), Float64(data["k_detrap"]))
        boundary = _build_boundary_operator_from_defs(bcs, build_hotgates_variable_layout(adaptor), [Float64(data["diffusion_coefficient"]), 0.0])
        return build_hotgates_trapping_model(
            mesh = mesh,
            model = backend,
            adaptor = adaptor,
            temperature = ConstantTemperature(Float64(data["temperature"])),
            diffusion_coefficients = [Float64(data["diffusion_coefficient"]), 0.0],
            boundary = boundary,
        )
    end

    return BackendDefinition(name, :hotgates_trapping, species, build_model, (;))
end

function _build_boundary_operator_from_defs(bcs::Vector{BCDefinition}, layout::VariableLayout, coefficients)
    isempty(bcs) && return nothing
    ops = AbstractOperator[]
    names = variable_names(layout)
    for bc in bcs
        idx = findfirst(==(bc.species), names)
        idx === nothing && throw(ConfigValidationError(
            "Boundary condition `$(bc.name)` references unknown species `$(bc.species)`. Available species: $(join(string.(names), ", "))",
        ))
        selector = let idx = idx
            layout -> [idx]
        end
        bc_fn = let value = bc.value
            t -> value
        end
        if bc.method == :weak
            if bc.boundary == :left
                push!(ops, WeakDirichletBoundaryOperator(selector, coefficients, nothing; left = bc_fn))
            else
                push!(ops, WeakDirichletBoundaryOperator(selector, coefficients, nothing; right = bc_fn))
            end
        else
            push!(ops, DirichletBoundaryOperator(
                bc.boundary,
                bc_fn,
                selector,
                coefficients,
                nothing;
                method = _boundary_method_from_symbol(bc.method),
            ))
        end
    end
    return length(ops) == 1 ? ops[1] : OperatorSum(Tuple(ops))
end

function _boundary_method_from_symbol(method::Symbol)
    if method == :penalty
        return PenaltyMethod()
    elseif method == :mass_matrix
        return MassMatrixMethod()
    elseif method == :callback
        return CallbackMethod()
    elseif method == :eliminated
        return EliminatedMethod()
    end
    throw(ArgumentError("Unsupported boundary method $method"))
end
