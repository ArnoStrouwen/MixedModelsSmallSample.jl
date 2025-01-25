using CSV
using DataFrames
using MixedModels
using LinearAlgebra
using ForwardDiff
using Distributions
using Test

using KenwardRoger

df = DataFrame(CSV.File("Data Pastry Dough Experiment Chapter 7.csv"))
rename!(df, "Flow Rate" => :FR, "Moisture Content" => :MC, "Screw Speed" => :SS)
rename!(df, "Longitudinal Expansion Index" => :LEI)
function one_minus_one_coding!(x)
    minx = minimum(x)
    maxx = maximum(x)
    span = maxx - minx
    x .-= minx .+ span / 2
    x ./= (span / 2)
    return nothing
end
one_minus_one_coding!(df[!, :FR])
one_minus_one_coding!(df[!, :MC])
one_minus_one_coding!(df[!, :SS])
fm = @formula(
    LEI ~
        1 +
        (1 | Day) +
        FR +
        MC +
        SS +
        FR & MC +
        FR & SS +
        MC & SS +
        FR & FR +
        MC & MC +
        SS & SS
)
m = fit(MixedModel, fm, df; REML=true)
kr = kenwardroger_matrices(m)
estimates = kenwardroger_estimates(m, kr)

res = DataFrame(CSV.File("Results pastry dough.csv"))
@test isapprox(res[!, "Estimate"], getfield.(estimates, :estimate), atol=1e-9, rtol=1e-9)
@test isapprox(res[!, "Std Error"], getfield.(estimates, :std_error), atol=1e-5, rtol=1e-6)
@test isapprox(res[!, "DFDen"], getfield.(estimates, :den_df), atol=1e-2, rtol=1e-4)
