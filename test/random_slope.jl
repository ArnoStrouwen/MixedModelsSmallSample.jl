using CSV
using DataFrames
using MixedModels
using Test

using KenwardRoger

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = kenwardroger_matrices(m)
estimates = kenwardroger_estimates(m, kr)

res = DataFrame(CSV.File("Results sleep study.csv"))
@test isapprox(res[!, "Estimate"], getfield.(estimates, :estimate), atol=1e-9, rtol=1e-9)
@test_broken isapprox(
    res[!, "Std Error"], getfield.(estimates, :std_error), atol=1e-5, rtol=1e-5
)
@test_broken isapprox(res[!, "DFDen"], getfield.(estimates, :den_df), atol=1e-5, rtol=1e-5)

kr = kenwardroger_matrices(m; FIM_σ²=:expected)
estimates = kenwardroger_estimates(m, kr)
res = DataFrame(CSV.File("Results sleep study lmertest.csv"))
@test isapprox(
    res[!, "coefficients.Estimate"], getfield.(estimates, :estimate), atol=1e-10, rtol=1e-10
)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    getfield.(estimates, :std_error),
    atol=1e-10,
    rtol=1e-8,
)
@test isapprox(
    res[!, "coefficients.df"], getfield.(estimates, :den_df), atol=1e-10, rtol=1e-9
)
