module Flopsy

include("config.jl")
include("types.jl")
include("mesh.jl")
include("state.jl")
include("reaction.jl")
include("diffusion.jl")
include("rhs.jl")
include("solver.jl")
include("output.jl")
include("runner.jl")

export run_simulation

end
