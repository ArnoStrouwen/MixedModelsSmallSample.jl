using Test
using DataFrames

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=false)
@test_throws AssertionError adjust_KR(m; FIM_σ²=:observed)
@test_throws AssertionError adjust_SW(m; FIM_σ²=:observed)
@test_throws AssertionError ftest_KR(m::LinearMixedModel, ones(2, 1))
@test_throws AssertionError ftest_SW(m::LinearMixedModel, ones(2, 1))
