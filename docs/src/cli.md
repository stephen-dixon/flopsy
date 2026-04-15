# CLI Usage

## Running From Source

The source-tree CLI wrapper is:

```bash
julia --project=. scripts/flopsy.jl --help
```

Main commands:

```text
flopsy run <input.toml>
flopsy validate <input.toml>
flopsy syntax list
flopsy syntax show <domain> <type>
flopsy plugin list
flopsy plugin register <name> [--registry <url>] [--url <pkg-url> | --path <local-path>]
flopsy plugin remove <name>
flopsy xdmf <fields.h5> [--output out.xdmf]
flopsy --version
```

## Main Workflows

Validate an input deck without running:

```bash
julia --project=. scripts/flopsy.jl validate examples/minimal_input_deck.toml
```

Run an input deck:

```bash
julia --project=. scripts/flopsy.jl run examples/minimal_input_deck.toml
```

Inspect syntax:

```bash
julia --project=. scripts/flopsy.jl syntax list
julia --project=. scripts/flopsy.jl syntax show problem simulation
```

Generate an XDMF companion for an existing HDF5 output:

```bash
julia --project=. scripts/flopsy.jl xdmf minimal_fields.h5
```

## Error Handling

By default the CLI reports user-facing validation and argument errors without a raw stack trace. For development/debugging, rerun with `--verbose` to include the full Julia error display.

Examples of improved validation coverage include:

- unknown keys within a block
- missing required fields
- invalid field types
- invalid enum-like values such as unsupported solver methods
- unknown syntax types
- invalid object references across named blocks
- invalid species targets for ICs and BCs

## Exit Codes

- `0` on success
- `1` on validation, argument, plugin-loading, or runtime errors

## Validate Before Run

If you want CI or scripted validation without solving:

```bash
julia --project=. scripts/flopsy.jl validate examples/trapping_1d_with_bc_ic.toml
```

This parses the deck, validates it against the live registry, builds the named objects, and assembles the referenced problem without entering the solver.
