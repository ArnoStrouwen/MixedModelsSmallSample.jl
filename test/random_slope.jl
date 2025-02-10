using CSV
using DataFrames
using MixedModels
using Test

using KenwardRoger

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = kenwardroger_matrices(m)
estimates = coeftable(m, kr)

res = DataFrame(CSV.File("Results sleep study.csv"))
@test isapprox(res[!, "Estimate"], estimates.cols[1], atol=1e-9, rtol=1e-9)
@test_broken isapprox(
    res[!, "Std Error"], estimates.cols[2], atol=1e-5, rtol=1e-5
)
@test_broken isapprox(res[!, "DFDen"], estimates.cols[6], atol=1e-5, rtol=1e-5)

kr = kenwardroger_matrices(m; FIM_σ²=:expected)
estimates = coeftable(m, kr)
res = DataFrame(CSV.File("Results sleep study lmertest.csv"))
@test isapprox(
    res[!, "coefficients.Estimate"], estimates.cols[1], atol=1e-10, rtol=1e-10
)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    estimates.cols[2],
    atol=1e-10,
    rtol=1e-8,
)
@test isapprox(
    res[!, "coefficients.df"], estimates.cols[6], atol=1e-10, rtol=1e-9
)
