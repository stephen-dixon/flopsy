using Flopsy

function main(args)
    length(args) == 1 ||
        error("usage: julia --project=. scripts/run_from_toml.jl <input.toml>")
    result = run_input_deck(args[1])
    println("retcode: ", getproperty(result.solution, :retcode))
    return result
end

main(ARGS)
