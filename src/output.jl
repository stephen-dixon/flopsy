module Output

using HDF5
using CSV
using DataFrames

export write_outputs

function write_outputs(sol, ctx, config)
    # --- CSV (observables placeholder)
    df = DataFrame(time=sol.t)
    CSV.write("observables.csv", df)

    # --- HDF5 full output
    h5open("fields.h5", "w") do f
        f["time"] = sol.t
        f["state"] = hcat(sol.u...)  # (state_dim, nt)
    end
end

end
