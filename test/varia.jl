using Test
using DataFrames

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels

df = DataFrame(MixedModels.dataset(:sleepstudy))
# MixedModelsDatasets >= 0.2 stores `reaction` as Float32; recover the exact
# 4-decimal values the reference results (SAS, JMP, lmerTest) were computed from.
df.reaction = parse.(Float64, string.(df.reaction))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=false)
@test_throws ErrorException small_sample_adjust(m, KenwardRoger(fim=ObservedFIM()))
small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))
@test_throws ErrorException small_sample_ftest(m, ones(2, 1); method=KenwardRoger())
small_sample_ftest(m, ones(2, 1); method=Satterthwaite())

@test_throws ErrorException small_sample_adjust(m, Satterthwaite(; fim=ExpectedFIM()))
@test_throws ErrorException vcov_varpar(m; fim=ExpectedFIM())
