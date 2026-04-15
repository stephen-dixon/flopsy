using Flopsy, Palioxis

Palioxis.init("/path/to/palioxis")
pal = Palioxis.MultipleDefectModel("model.xml")

mesh = Mesh1D(1e-3, 200)
nx = length(mesh.x)

# Build per-node defect profile from the model's XML depth profile
defects = hcat([Palioxis.get_defect_concentrations(pal, mesh.x[i] - mesh.dx/2, mesh.x[i] +
                                                                               mesh.dx/2)
                for i in 1:nx]...)

# D(T) is queried from Palioxis at every time step automatically
model = build_palioxis_trapping_model(
    palioxis_model = pal,
    mesh = mesh,
    defects = defects,
    temperature = ConstantTemperature(600.0)
)

# With surface boundary conditions:
# model = build_palioxis_trapping_model(
#     palioxis_model = pal,
#     mesh           = mesh,
#     defects        = defects,
#     temperature    = LinearRampTemperature(300.0, 10/60),
#     left_bc        = t -> 0.0,
#     right_bc       = t -> 0.0,
# )
