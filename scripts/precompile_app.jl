using Flopsy

Flopsy.cli_main(["--help"])
Flopsy.cli_main(["syntax", "list"])
Flopsy.cli_main(["validate", joinpath(@__DIR__, "..", "examples", "minimal_input_deck.toml")])
