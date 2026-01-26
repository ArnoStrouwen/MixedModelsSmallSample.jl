using Random
using LinearAlgebra

using CSV
using DataFrames

function generate_fa0_drop_simulation(;
    seed::Integer=1,
    rhos::Vector{Float64}=[
        0.90, 0.95, 0.97, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9995, 0.9999, 0.99995, 0.99999
    ],
    nrep::Integer=50,
    nsubj::Integer=24,
    times::Vector{Float64}=[-1.0, 0.0, 1.0],
    beta0::Float64=0.0,
    beta1::Float64=1.0,
    sd_intercept::Float64=1.0,
    sd_slope::Float64=0.05,
    sd_eps::Float64=0.01,
    out_csv::AbstractString=joinpath(
        @__DIR__, "..", "data", "Data fa0 drop simulation.csv"
    ),
)
    rng = MersenneTwister(seed)

    rows = NamedTuple[]
    sim_id = 0

    for rho in rhos
        G = [
            sd_intercept^2 rho * sd_intercept * sd_slope
            rho * sd_intercept * sd_slope sd_slope^2
        ]
        L = cholesky(Symmetric(G); check=false).L

        for rep in 1:nrep
            sim_id += 1
            for subj in 1:nsubj
                z = randn(rng, 2)
                b0, b1 = L * z

                for t in times
                    y = (beta0 + beta1 * t) + (b0 + b1 * t) + sd_eps * randn(rng)
                    push!(
                        rows, (sim_id=sim_id, rho=rho, rep=rep, subject=subj, time=t, y=y)
                    )
                end
            end
        end
    end

    df = DataFrame(rows)
    CSV.write(out_csv, df)
    return out_csv
end

if abspath(PROGRAM_FILE) == @__FILE__
    out = generate_fa0_drop_simulation()
    println("Wrote: ", out)
end
