using Flopsy

deck_path = joinpath(@__DIR__, "minimal_input_deck.toml")

ctx = validate_input_deck(deck_path)
println("validated problems: ", join(sort!(string.(collect(keys(ctx.problems)))), ", "))

result = run_input_deck(deck_path)
println("retcode: ", getproperty(result.solution, :retcode))
