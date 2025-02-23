using CSV
using DataFrames
using MixedModels
using Test
using LinearAlgebra

using KenwardRoger

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = adjust_KR(m; FIM_σ²=:observed_SAS_MATCHING)

res = DataFrame(CSV.File("Results sleep study jmp.csv"))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(
    res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-6
)
@test isapprox(res[!, "DFDen"], kr.v, atol=1e-10, rtol=1e-5)

res = DataFrame(CSV.File("Results sleep study sas kr.csv"))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-5, rtol=1e-10)
@test isapprox(res[!, "DF"], kr.v, atol=1e-10, rtol=1e-5)

sw = adjust_SW(m; FIM_σ²=:observed_SAS_MATCHING)

res = DataFrame(CSV.File("Results sleep study sas sw.csv"))
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-5, rtol=1e-10)
@test isapprox(res[!, "DF"], sw.v, atol=1e-10, rtol=1e-5)

kr = adjust_KR(m; FIM_σ²=:expected)

res = DataFrame(CSV.File("Results sleep study lmertest.csv"))
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-10,
    rtol=1e-8,
)
@test isapprox(res[!, "coefficients.df"], kr.v, atol=1e-10, rtol=1e-9)
