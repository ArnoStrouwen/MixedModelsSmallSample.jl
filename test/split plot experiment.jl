using CSV
using DataFrames
using MixedModels
using Test
using LinearAlgebra

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

kr = adjust_KR(m; FIM_σ²=:observed_SAS_MATCHING)

res = DataFrame(CSV.File("Results wind tunnel.csv"))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-9, rtol=1e-9)
@test isapprox(
    res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-7
)
@test isapprox(res[!, "DFDen"], kr.v, atol=1e-10, rtol=1e-8)

kr = adjust_KR(m; FIM_σ²=:expected)

res = DataFrame(CSV.File("Results wind tunnel lmertest.csv"))
res = vcat(res, res[6:9, :])
deleteat!(res, 6:9)
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-6,
    rtol=1e-10,
)
@test isapprox(res[!, "coefficients.df"], kr.v, atol=1e-9, rtol=1e-9)
