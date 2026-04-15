# Flopsy.jl Documentation

## What Flopsy Is

Flopsy is a Julia package for stiff one-dimensional reaction-diffusion simulations.

It has two supported entry styles:

- a CLI for registry-driven TOML input decks
- a Julia library API for direct model construction and execution

These are not separate numerical engines. The input-deck system validates and assembles objects that ultimately run through the same solver stack used by the Julia API.

## Supported TOML Workflow

The only supported TOML configuration workflow is the registry-driven input-deck model.

Input decks are composed from named blocks in fixed domains:

- `mesh`
- `backend`
- `ic`
- `bc`
- `output`
- `problem`

Each block declares a registered `type`. Built-ins and plugins both contribute syntax entries through the same `SyntaxRegistry`.

The old typed TOML/config path still exists only as a deprecated compatibility layer and is no longer documented as a primary workflow.

## Primary User Paths

### CLI

Use the CLI when you want:

- a declarative input-deck workflow
- centralized validation and diagnostics
- live syntax inspection
- runtime plugin loading
- a clean path to a standalone app build

### Julia Library

Use the Julia API when you want:

- direct access to `SystemModel` construction
- programmatic composition of simulations
- custom orchestration from Julia code
- lower-level extension points

## Canonical Flow

The canonical input-deck flow is:

1. Parse the input deck.
2. Validate blocks against registered syntax.
3. Build a validated `BuildContext` of named objects.
4. Assemble a `ConfiguredSimulation` from a `[problem.<name>]` block.
5. Run it through the CLI or from Julia.

## Quick Links

- [CLI Usage](cli.md)
- [Julia Library Usage](julia_library.md)
- [Input Deck Format](input_deck.md)
- [Syntax Registration and Extension Model](extensions.md)
- [Plugins](plugins.md)
- [Standalone Compilation](standalone_compilation.md)
- [Legacy API Deprecation Note](legacy_api.md)
- [Architecture](architecture.md)
- [API Reference](api.md)
