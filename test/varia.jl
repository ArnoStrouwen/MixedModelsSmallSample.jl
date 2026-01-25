using Test
using DataFrames

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=false)
@test_throws ErrorException small_sample_adjust(m, KenwardRoger(fim=ObservedFIM()))
small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))
@test_throws ErrorException small_sample_ftest(m, ones(2, 1); method=KenwardRoger())
small_sample_ftest(m, ones(2, 1); method=Satterthwaite())

@test_throws ErrorException small_sample_adjust(m, Satterthwaite(; fim=ExpectedFIM()))
@test_throws ErrorException vcov_varpar(m; fim=ExpectedFIM())
