using Documenter
using Flopsy

makedocs(
    sitename = "Flopsy.jl",
    format = Documenter.HTML(edit_link = "main"),
    modules = [Flopsy],
    pages = [
        "Overview" => "index.md",
        "CLI Usage" => "cli.md",
        "Julia Library Usage" => "julia_library.md",
        "Input Deck Format" => "input_deck.md",
        "Syntax Registration" => "extensions.md",
        "Plugins" => "plugins.md",
        "Standalone Compilation" => "standalone_compilation.md",
        "Legacy API" => "legacy_api.md",
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
