using Flopsy

function main(args)
    length(args) == 1 || error("usage: julia --project=. scripts/run_from_toml.jl <input.toml>")
    cfg = load_config(args[1])
    problem = build_problem(cfg)
    result = solve(problem)
    println("retcode: ", getproperty(result.solution, :retcode))
    return result
end

main(ARGS)
