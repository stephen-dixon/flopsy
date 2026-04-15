# Plugins

## Plugin Contract

Flopsy plugins extend the input-deck syntax registry.

A runtime plugin package must provide:

```julia
register_flopsy_plugin!(registry)
```

Inside that function, the plugin registers additional `SyntaxSpec` entries.

## Runtime Loading Model

The CLI discovers runtime plugins from a managed Julia environment under the user config directory.

This keeps the compiled CLI stable while still allowing new syntax to be installed after the app is built.

Operationally:

1. Flopsy creates or reuses the managed plugin environment.
2. Registered plugin packages are loaded from that environment.
3. Each plugin's `register_flopsy_plugin!(registry)` function is invoked.
4. Plugin load failures are recorded and surfaced through `plugin list` and CLI errors.

## User Workflow

List installed plugins:

```bash
julia --project=. scripts/flopsy.jl plugin list
```

Register a local plugin:

```bash
julia --project=. scripts/flopsy.jl plugin register MyPlugin --path /path/to/MyPlugin
```

Register from a package URL:

```bash
julia --project=. scripts/flopsy.jl plugin register MyPlugin --url https://example.com/MyPlugin.jl.git
```

Remove a plugin:

```bash
julia --project=. scripts/flopsy.jl plugin remove MyPlugin
```

## What Users Write in TOML

Once a plugin is installed and loaded successfully, its registered syntax appears alongside built-ins in:

- `syntax list`
- `syntax show <domain> <type>`

Users then reference the plugin's registered `(domain, type)` pairs in normal input-deck blocks. Plugins do not invent new top-level TOML domains.

## Built-In Extensions

The Palioxis integration in `ext/PalioxisExt.jl` follows the same registration model conceptually:

- it contributes backend/IC syntax through the registry
- it coexists with built-in syntax
- if unavailable, its syntax is simply absent rather than partially parsed by a special TOML path
