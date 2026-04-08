using Documenter
using Flopsy

makedocs(
    sitename = "Flopsy.jl",
    format   = Documenter.HTML(),
    modules  = [Flopsy],
    pages    = [
        "Home"              => "index.md",
        "Architecture"      => "architecture.md",
        "Hotgates Interface" => "hotgates.md",
        "API Reference"     => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/palioxis-tds/flopsy.git",
    devbranch = "main",
)
