# Julia Library Usage

## Overview

Flopsy is not only a CLI tool. It is also a Julia library for constructing and running simulations directly.

There are two common Julia-side workflows:

- run a registry-driven input deck from Julia
- build the model and solver inputs programmatically

## Running an Input Deck From Julia

```julia
using Flopsy

deck_path = "examples/minimal_input_deck.toml"

ctx = validate_input_deck(deck_path)
result = run_input_deck(deck_path)
```

This gives you the same validation and assembly path that the CLI uses, but keeps control inside your Julia process.

If you want the intermediate build products:

```julia
using Flopsy

deck = parse_input_deck("examples/minimal_input_deck.toml")
registry = build_registry()
ctx = build_context(deck; registry = registry)
plan = build_simulation(ctx)
result = solve(plan)
```

## Building a Simulation Directly

```julia
using Flopsy
using OrdinaryDiffEq

mesh = Mesh1D(1.0, 51)
model = build_trapping_model(
    mesh = mesh,
    k_trap = 2.0,
    k_detrap = 0.25,
    diffusion_coefficient = 0.05,
)

u0 = zeros(nvariables(model.layout) * length(mesh.x))
tspan = (0.0, 0.25)

solver = SolverConfig(
    formulation = UnsplitFormulation(),
    algorithm = Rodas5P(),
    abstol = 1e-8,
    reltol = 1e-6,
    saveat = [0.0, 0.125, 0.25],
)

result = wrap_result(model, solve_problem(model, u0, tspan, solver), solver)
```

This is the right path when you want full programmatic control over model construction, solver choice, orchestration, or post-processing.

## Lower-Level Extension Points

Advanced users can work below the input-deck layer by extending or composing:

- `SystemModel`
- `VariableLayout`
- reaction, diffusion, and constraint operators
- formulations and solver configuration
- syntax registration through `SyntaxRegistry`

## When To Prefer Each Path

Prefer the CLI or `run_input_deck` when:

- you want a reproducible declarative deck
- you want centralized schema validation
- you want plugin-discovered syntax

Prefer direct Julia construction when:

- your setup is naturally programmatic
- you are integrating Flopsy into a larger Julia workflow
- you need lower-level control of model objects or solve orchestration
