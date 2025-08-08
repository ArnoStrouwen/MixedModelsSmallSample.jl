using MixedModels
using Test
using DataFrames

using MixedModelsSmallSample

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=false)
@test_throws AssertionError adjust(m; method=KenwardRoger(), FIM_σ²=:observed_SAS_MATCHING)
@test_throws AssertionError adjust(m; method=Satterthwaite(), FIM_σ²=:observed_SAS_MATCHING)
@test_throws AssertionError ftest_KR(m::LinearMixedModel, ones(2, 1))
@test_throws AssertionError ftest_SW(m::LinearMixedModel, ones(2, 1))

# Test new adjust interface with REML=false (should also throw AssertionError)
@test_throws AssertionError adjust(m; method=KenwardRoger())
@test_throws AssertionError adjust(m; method=Satterthwaite())
@test_throws AssertionError adjust(m)  # Test default method
