using Flopsy, Palioxis

  Palioxis.init("/path/to/palioxis")
  pal = Palioxis.MultipleDefectModel("model.xml")

  # build per-node defect profile from the model's XML depth profile
  defects = hcat([Palioxis.get_defect_concentrations(pal, mesh.x[i]-mesh.dx/2, mesh.x[i]+mesh.dx/2) for i in 1:nx]...)

  model = build_palioxis_trapping_model(
      palioxis_model        = pal,
      mesh                  = mesh,
      defects               = defects,
      temperature           = ConstantTemperature(600.0),
      diffusion_coefficients = Palioxis.diffusion_constants(pal, 600.0),
  )
