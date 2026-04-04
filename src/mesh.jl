function Mesh1D(length::Real, nx::Integer)
    nx >= 2 || throw(ArgumentError("nx must be at least 2"))
    x = collect(range(zero(length), length; length=nx))
    dx = x[2] - x[1]
    return Mesh1D(x, dx)
end
