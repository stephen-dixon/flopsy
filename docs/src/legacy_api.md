# Legacy Configuration API

The old typed config/TOML path is deprecated and is no longer the recommended user workflow.

Deprecated compatibility entry points include:

- `load_config`
- `parse_config`
- `validate(::ProblemConfig)`
- `build_problem(::ProblemConfig)`
- the old typed config structs such as `ProblemConfig`

## Current Status

- the code remains temporarily for compatibility
- deprecation warnings are emitted when those entry points are used
- README and primary docs no longer present that path as the main TOML workflow
- new TOML-facing features should be added only through the registry-driven input-deck system

## Supported TOML Path

The supported TOML configuration path is:

1. `parse_input_deck`
2. `build_context`
3. `build_simulation`
4. `solve` or `run_input_deck`

If you still depend on the typed config layer, treat it as a migration target rather than a stable forward-looking API.
