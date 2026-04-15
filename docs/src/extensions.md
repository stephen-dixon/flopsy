# Syntax Registration and Extension Model

## Registry Overview

The registry is the canonical TOML extension point.

Each syntax entry is registered by `(domain, type)` and carries:

- parameter schema (`ParameterSpec`)
- defaults
- required/optional status
- type expectations
- optional allowed values
- help text
- validation callback
- build callback
- plugin/source attribution

Built-in syntax and plugin syntax use the same mechanism.

## ParameterSpec

`ParameterSpec` is the central schema object used for user-facing TOML validation.

It can describe:

- field name
- whether the field is required
- default value
- field kind such as `:string`, `:real`, `:integer`, `:bool`, or `:vector`
- allowed enum-like values
- doc/help text

That schema is enforced centrally during `build_context(deck)` before builders execute.

## Registration Flow

A syntax provider registers entries into a `SyntaxRegistry`:

```julia
register_syntax!(registry, SyntaxSpec(
    :backend,
    :my_backend,
    [
        ParameterSpec(:type, true, nothing, "Syntax type"; kind = :string),
        ParameterSpec(:coefficient, true, nothing, "Transport coefficient"; kind = :real),
    ],
    "Custom backend example.",
    validate_fn,
    build_fn,
    :my_plugin,
))
```

## Builder Responsibilities

The schema layer should catch structural TOML issues before builder code runs.

Builder code should focus on:

- domain-specific semantic validation not expressible in the generic schema
- constructing the runtime objects
- reporting unresolved references or unsupported targets clearly

## Built-In and Plugin Syntax Together

At runtime Flopsy:

1. registers built-in syntax
2. loads runtime plugins from the managed plugin environment
3. lets each plugin register additional syntax into the same registry

As a result, `syntax list` and `syntax show` expose one combined view.

## Experimental or Unavailable Syntax

If a syntax family is registered before its implementation is production-ready, it should be clearly marked as experimental and should fail explicitly rather than silently appearing supported.

In the current tree, `biased_1d` is discoverable for architectural visibility but reports that it is experimental and not ready for production use.
