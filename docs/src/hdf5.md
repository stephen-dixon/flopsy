# HDF5 Output

Flopsy writes field output to HDF5 files and can generate an XDMF companion so the results can be opened directly in ParaView.

## Producing output from TOML

The built-in output syntax is:

```toml
[output.fields]
type = "hdf5"
file = "fields.h5"
xdmf = true
```

If `xdmf = true`, Flopsy writes the HDF5 file and then emits a matching `.xdmf` companion automatically.

You can also generate the companion afterwards:

```bash
flopsy xdmf fields.h5
flopsy xdmf fields.h5 --output custom_name.xdmf
```

## ParaView workflow

For ParaView:

1. run the simulation and produce `fields.h5`
2. generate or request the matching `.xdmf`
3. open the `.xdmf` file in ParaView

The XDMF file references the HDF5 datasets rather than duplicating the field data.

## File Structure

Files written by `write_field_output_hdf5(result, path)` follow this layout:

```
/time              — Float64 vector of saved simulation times (length nt)
/mesh/x            — Float64 vector of spatial node coordinates (length nx)
/mesh/dx           — Float64 scalar, uniform node spacing
/fields/<varname>  — Float64 matrix of shape (nt, nx), one entry per variable
                      e.g. mobile_H, trap_H_1, trap_H_2, ...
/metadata/...      — Selected simulation metadata stored as string datasets
                      e.g. model_type, xml_file, retcode, walltime_s, ...
```

Variable names are Julia Symbols stored as plain strings in the HDF5 file. The field matrices are indexed `[time_index, node_index]`.

## Reading in Julia

Use [HDF5.jl](https://github.com/JuliaIO/HDF5.jl):

```julia
using HDF5

h5open("tds_fields.h5", "r") do h5
    t  = read(h5["time"])          # Vector{Float64}, length nt
    x  = read(h5["mesh/x"])        # Vector{Float64}, length nx
    dx = read(h5["mesh/dx"])       # Float64

    # List available field variables
    vars = keys(h5["fields"])      # e.g. ["mobile_H", "trap_H_1", ...]

    # Read one variable — shape (nt, nx)
    mobile_H = read(h5["fields"]["mobile_H"])

    # Spatial profile at the final saved time
    profile_end = mobile_H[end, :]

    # Spatial integral over time (trapezoidal with uniform dx)
    inventory = sum(mobile_H, dims=2)[:] .* dx
end
```

### Read all fields into a Dict

```julia
using HDF5

h5open("tds_fields.h5", "r") do h5
    t  = read(h5["time"])
    x  = read(h5["mesh/x"])
    dx = read(h5["mesh/dx"])

    fields = Dict(varname => read(h5["fields"][varname])
                  for varname in keys(h5["fields"]))

    metadata = Dict(k => read(h5["metadata"][k])
                    for k in keys(h5["metadata"]))

    for (name, data) in fields
        println("$name  shape=$(size(data))")
    end
end
```

### Desorption flux from a TDS run

The HDF5 file stores field profiles; surface fluxes are most easily obtained from the
summary CSV. If you need them from the HDF5 directly:

```julia
using HDF5

h5open("tds_fields.h5", "r") do h5
    t        = read(h5["time"])
    dx       = read(h5["mesh/dx"])
    # Example: constant D for illustration — replace with actual D(T) if needed
    D        = 1e-7
    mobile_H = read(h5["fields"]["mobile_H"])   # (nt, nx)

    # Right-surface flux (outward normal, positive during desorption)
    nx = size(mobile_H, 2)
    right_flux = D .* (mobile_H[:, nx-1] .- mobile_H[:, nx]) ./ dx
end
```

### Chaining simulations

`load_ic_from_hdf5` reads the HDF5 output from one run and constructs an initial-condition
vector for a subsequent simulation.  Variables are matched by name; missing variables are
zero-initialised with a warning.

```julia
using Flopsy

# --- Stage 1: implantation ---
sol_impl = solve_problem(model_impl, u0_impl, tspan_impl, config_impl)
result_impl = wrap_result(model_impl, sol_impl, config_impl)
write_field_output_hdf5(result_impl, "implant.h5")

# --- Stage 2: TDS ramp starting from the implanted state ---
u0_tds = load_ic_from_hdf5("implant.h5", model_tds)               # load final state
# or: load_ic_from_hdf5("implant.h5", model_tds; time_index = 10) # load nth saved time

sol_tds = solve_problem(model_tds, u0_tds, tspan_tds, config_tds)
result_tds = wrap_result(model_tds, sol_tds, config_tds)
```

## Reading in Python

Use [h5py](https://www.h5py.org/) and NumPy:

```python
import h5py
import numpy as np

with h5py.File("tds_fields.h5", "r") as h5:
    t  = h5["time"][:]          # 1-D array, shape (nt,)
    x  = h5["mesh/x"][:]        # 1-D array, shape (nx,)
    dx = float(h5["mesh/dx"][()])

    # List available variables
    var_names = list(h5["fields"].keys())
    print("Variables:", var_names)

    # Read one variable — shape (nt, nx)
    mobile_H = h5["fields"]["mobile_H"][:]

    # Final profile and total inventory
    profile_end = mobile_H[-1, :]                    # last time step
    inventory   = mobile_H.sum(axis=1) * dx          # scalar per time step

    # Read string metadata
    model_type = h5["metadata"]["model_type"][()].decode()
    retcode    = h5["metadata"]["retcode"][()].decode()
```

!!! note "Index ordering"
    Julia arrays are 1-indexed and column-major.  HDF5.jl writes `(nt, nx)` matrices in
    row-major order on disk (time is the slow axis) to match the convention used by h5py
    and NumPy.  Python code reading `h5["fields"]["mobile_H"][:]` will always get shape
    `(nt, nx)` with `array[-1, :]` being the last time step.

### Plot desorption flux vs temperature (Python)

The summary CSV is the easiest source for flux vs T plots.  With the HDF5 file directly:

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Suppose temperature was saved as metadata in a companion CSV or extra field.
# Here we reconstruct a linear ramp as an example.

with h5py.File("tds_fields.h5", "r") as h5:
    t        = h5["time"][:]
    dx       = float(h5["mesh/dx"][()])
    mobile_H = h5["fields"]["mobile_H"][:]

nx = mobile_H.shape[1]

# Constant D for illustration
D = 1e-7
right_flux = D * (mobile_H[:, -2] - mobile_H[:, -1]) / dx

# Linear temperature ramp: 300 K + 10 K/min
T_ramp = 300.0 + (10.0 / 60.0) * t

plt.figure(figsize=(7, 4))
plt.plot(T_ramp, right_flux, lw=1.5)
plt.xlabel("Temperature (K)")
plt.ylabel("Right-surface flux (m$^{-2}$ s$^{-1}$)")
plt.title("TDS desorption spectrum")
plt.tight_layout()
plt.savefig("tds_spectrum.png", dpi=150)
```

## XDMF helper

Flopsy provides a native helper:

```julia
using Flopsy

write_xdmf_for_flopsy_h5("fields.h5")
write_xdmf_for_flopsy_h5("fields.h5", "fields_custom.xdmf")
```

The CLI `flopsy xdmf` subcommand is a thin wrapper over the same functionality.

## Output functions reference

| Function | Output | Use |
|---|---|---|
| `write_field_output_hdf5(result, path)` | HDF5 file with full `(nt, nx)` fields | Spatial analysis, chaining |
| `write_summary_csv(result, path)` | CSV with integrals and surface fluxes per time step | Quick inspection, flux plots |
| `load_ic_from_hdf5(path, model; time_index=:last)` | Flat `Vector{Float64}` state | Chaining simulations |

### Summary CSV columns

The CSV produced by `write_summary_csv` contains:

| Column | Description |
|---|---|
| `time` | Saved simulation time |
| `integral_<varname>` | Spatial integral of each variable (one column per variable) |
| `left_flux_<varname>` | Outward diffusive flux at left boundary (x = 0) |
| `right_flux_<varname>` | Outward diffusive flux at right boundary (x = L) |
| extra columns | Any entries in `result.summaries[:extra_timeseries]`, e.g. `temperature_K` |

Both fluxes are positive when hydrogen leaves the sample (vacuum boundary condition at TDS surfaces).
