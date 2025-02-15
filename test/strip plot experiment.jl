using CSV
using DataFrames
using MixedModels
using Test

using KenwardRoger

df = DataFrame(CSV.File("Data battery cell Chapter 11.csv"))
rename!(df, "Whole Plots" => :WP)
rename!(df, "Subplots" => :SP)

fm = @formula(
    Y ~
        1 +
        (1 | WP) +
        (1 | SP) +
        X1 +
        X2 +
        X3 +
        X4 +
        X5 +
        X6 +
        (X2 + X3 + X4 + X5 + X6) & (X1) +
        (X3 + X4 + X5 + X6) & (X2) +
        (X4 + X5 + X6) & (X3) +
        (X5 + X6) & (X4) +
        (X6) & (X5)
)
m = fit(MixedModel, fm, df; REML=true)

kr = kenwardroger_matrices(m)
estimates = coeftable(m, kr)

res = DataFrame(CSV.File("Results battery cell.csv"))
@test isapprox(res[!, "Estimate"], estimates.cols[1], atol=1e-7, rtol=1e-7)
@test_broken isapprox(res[!, "Std Error"], estimates.cols[2], atol=1e-5, rtol=1e-5)
@test_broken isapprox(res[!, "DFDen"], estimates.cols[6], atol=1e-5, rtol=1e-5)

kr = kenwardroger_matrices(m; FIM_σ²=:expected)
estimates = coeftable(m, kr)

res = DataFrame(CSV.File("Results battery cell lmertest.csv"))
@test isapprox(res[!, "coefficients.Estimate"], estimates.cols[1], atol=1e-10, rtol=1e-7)
@test isapprox(res[!, "coefficients.Std..Error"], estimates.cols[2], atol=1e-5, rtol=1e-10)
@test isapprox(res[!, "coefficients.df"], estimates.cols[6], atol=1e-10, rtol=1e-6)
