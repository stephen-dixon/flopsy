using Documenter
using Flopsy

makedocs(
    sitename = "Flopsy.jl",
    format = Documenter.HTML(edit_link = "main"),
    modules = [Flopsy],
    pages = [
        "Home" => "index.md",
        "Configuration" => "configuration.md",
        "Architecture" => "architecture.md",
        "Formulations" => "formulations.md",
        "Hotgates Interface" => "hotgates.md",
        "HDF5 Output" => "hdf5.md",
        "Plotting" => "plotting.md",
        "API Reference" => "api.md"
    ],
    checkdocs = :exports
)

deploydocs(
    repo = "github.com/stephen-dixon/flopsy.git",
    devbranch = "main"
)
