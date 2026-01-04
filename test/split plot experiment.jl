using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(CSV.File("data/Data wind tunnel Chapter 10.csv"))
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

kr = adjust_KR(m; FIM_σ²=:observed_SAS_MATCHING)

res = DataFrame(CSV.File("results/Results wind tunnel jmp.csv"))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-9, rtol=1e-8)
@test isapprox(
    res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-5
)
@test isapprox(res[!, "DFDen"], kr.v, atol=1e-10, rtol=1e-6)

res = DataFrame(CSV.File("results/Results wind tunnel sas kr.csv"))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-9, rtol=1e-8)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-7, rtol=1e-10)
@test isapprox(res[!, "DF"], kr.v, atol=1e-10, rtol=1e-6)

sw = adjust_SW(m; FIM_σ²=:observed_SAS_MATCHING)

res = DataFrame(CSV.File("results/Results wind tunnel sas sw.csv"))
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-9, rtol=1e-8)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-7, rtol=1e-10)
@test isapprox(res[!, "DF"], sw.v, atol=1e-10, rtol=1e-6)

kr = adjust_KR(m; FIM_σ²=:expected)

res = DataFrame(CSV.File("results/Results wind tunnel lmertest.csv"))
res = vcat(res, res[6:9, :])
deleteat!(res, 6:9)
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-10, rtol=1e-8)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-6,
    rtol=1e-9,
)
@test isapprox(res[!, "coefficients.df"], kr.v, atol=1e-9, rtol=1e-7)
