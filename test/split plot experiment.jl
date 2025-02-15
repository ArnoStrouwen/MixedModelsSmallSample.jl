using CSV
using DataFrames
using MixedModels
using Test

using KenwardRoger

df = DataFrame(CSV.File("Data wind tunnel Chapter 10.csv"))
rename!(df, "Whole Plots" => :WP)

fm = @formula(
    EFFICIENCY ~
        1 +
        (1 | WP) +
        FRH +
        RRH +
        YA +
        GC +
        FRH & RRH +
        FRH & YA +
        FRH & GC +
        RRH & YA +
        RRH & GC +
        YA & GC +
        FRH & FRH +
        RRH & RRH +
        YA & YA +
        GC & GC
)
m = fit(MixedModel, fm, df; REML=true)

kr = kenwardroger_matrices(m)
estimates = coeftable(m, kr)

res = DataFrame(CSV.File("Results wind tunnel.csv"))
@test isapprox(res[!, "Estimate"], estimates.cols[1], atol=1e-9, rtol=1e-9)
@test isapprox(res[!, "Std Error"], estimates.cols[2], atol=1e-6, rtol=1e-6)
@test isapprox(res[!, "DFDen"], estimates.cols[6], atol=1e-2, rtol=1e-2)

kr = kenwardroger_matrices(m; FIM_σ²=:expected)
estimates = coeftable(m, kr)

res = DataFrame(CSV.File("Results wind tunnel lmertest.csv"))
res = vcat(res, res[6:9, :])
deleteat!(res, 6:9)
@test isapprox(res[!, "coefficients.Estimate"], estimates.cols[1], atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "coefficients.Std..Error"], estimates.cols[2], atol=1e-6, rtol=1e-10)
@test isapprox(res[!, "coefficients.df"], estimates.cols[6], atol=1e-9, rtol=1e-9)
