using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
const APP_DIR = joinpath(ROOT, "app", "FlopsyCLI")
const BUILD_DIR = joinpath(ROOT, "build", "flopsy-app")

Pkg.activate(APP_DIR; io = devnull)
Pkg.develop(Pkg.PackageSpec(path = ROOT); io = devnull)

using PackageCompiler

create_app(
    APP_DIR,
    BUILD_DIR;
    executables = ["flopsy" => "julia_main"],
    precompile_execution_file = joinpath(ROOT, "scripts", "precompile_app.jl"),
    force = true
)

println("Built app at ", BUILD_DIR)
