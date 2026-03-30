module Mesh

export create_mesh

function create_mesh(nx::Int, length::Float64)
    x = range(0, length; length=nx)
    dx = step(x)
    return x, dx
end

end
