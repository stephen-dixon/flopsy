# Standalone Compilation

## Goal

Flopsy is structured so the CLI works:

- directly from source
- as a PackageCompiler app with a `julia_main()::Cint` entry point

## Entry Points

Source CLI wrapper:

```bash
julia --project=. scripts/flopsy.jl --help
```

App-facing entry point:

- package: `app/FlopsyCLI`
- function: `FlopsyCLI.julia_main()::Cint`

The wrapper delegates to `Flopsy.julia_main()`, which in turn calls the same CLI dispatch path used in normal Julia execution.

## Build Script

The repository includes:

- `scripts/build_app.jl`
- `scripts/precompile_app.jl`

Build the app with:

```bash
julia --project=. scripts/build_app.jl
```

Expected output:

- directory: `build/flopsy-app`
- executable name: `flopsy`

## PackageCompiler Notes

The app build is designed to avoid assumptions that only hold in an interactive development session:

- the CLI has a clean `julia_main()::Cint` entry point
- source paths are not embedded in user documentation or examples
- plugin loading is routed through a managed runtime environment instead of requiring the compiled app to rebuild itself

Relocatability still depends on downstream dependencies behaving well with artifacts/preferences. If additional packages start embedding absolute install-time paths into compiled code, that should be treated as a packaging bug and fixed at the dependency boundary.
