using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(CSV.File("data/Data bioequivalence.csv"))

fm = @formula(log_data ~ 1 + (1 | subject) + formulation + sequence + period)
m = fit(MixedModel, fm, df; REML=true, contrasts=Dict(:period => DummyCoding()))

res = DataFrame(CSV.File("results/Results bioequivalence jmp.csv"))

q = 1
L = zeros(length(m.betas), q)
L[2, 1] = 1
F_test_formulation = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[1, "DFDen"], F_test_formulation[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[1, "F Ratio"], F_test_formulation[2], atol=1e-10, rtol=1e-4)

q = 1
L = zeros(length(m.betas), q)
L[3, 1] = 1
F_test_sequence = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[3, "DFDen"], F_test_sequence[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[3, "F Ratio"], F_test_sequence[2], atol=1e-10, rtol=1e-5)

q = 3
L = zeros(length(m.betas), q)
L[4, 1] = 1
L[5, 2] = 1
L[6, 3] = 1
F_test_period = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[2, "DFDen"], F_test_period[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[2, "F Ratio"], F_test_period[2], atol=1e-10, rtol=1e-4)

res = DataFrame(CSV.File("results/Results bioequivalence homogeneous sas kr.csv"))

q = 1
L = zeros(length(m.betas), q)
L[2, 1] = 1
F_test_formulation = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[2, "DenDF"], F_test_formulation[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[2, "FValue"], F_test_formulation[2], atol=1e-10, rtol=1e-4)

q = 1
L = zeros(length(m.betas), q)
L[3, 1] = 1
F_test_sequence = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[3, "DenDF"], F_test_sequence[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[3, "FValue"], F_test_sequence[2], atol=1e-10, rtol=1e-5)

q = 3
L = zeros(length(m.betas), q)
L[4, 1] = 1
L[5, 2] = 1
L[6, 3] = 1
F_test_period = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[1, "DenDF"], F_test_period[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[1, "FValue"], F_test_period[2], atol=1e-10, rtol=1e-4)

res = DataFrame(CSV.File("results/Results bioequivalence homogeneous sas sw.csv"))

q = 1
L = zeros(length(m.betas), q)
L[2, 1] = 1
F_test_formulation = ftest_SW(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[2, "DenDF"], F_test_formulation[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[2, "FValue"], F_test_formulation[2], atol=1e-10, rtol=1e-6)

q = 1
L = zeros(length(m.betas), q)
L[3, 1] = 1
F_test_sequence = ftest_SW(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[3, "DenDF"], F_test_sequence[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[3, "FValue"], F_test_sequence[2], atol=1e-10, rtol=1e-6)

q = 3
L = zeros(length(m.betas), q)
L[4, 1] = 1
L[5, 2] = 1
L[6, 3] = 1
F_test_period = ftest_SW(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[1, "DenDF"], F_test_period[1], atol=1e-10, rtol=1e-5)
@test isapprox(res[1, "FValue"], F_test_period[2], atol=1e-10, rtol=1e-6)

#= df = transform(df, :log_data => eachindex => :row)

fm =  @formula(log_data ~ 1 + formulation + sequence + period + (formulation + 0 | subject) + zerocorr(formulation + 0 | row))
m = fit(MixedModel, fm, df; REML=true, contrasts=Dict(:period => DummyCoding()))

res = DataFrame(CSV.File("results/Results bioequivalence heterogeneous sas kr.csv"))

q = 1
L = zeros(length(m.betas), q)
L[2, 1] = 1
F_test_formulation = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[2, "DenDF"], F_test_formulation[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[2, "FValue"], F_test_formulation[2], atol=1e-10, rtol=1e-4)

q = 1
L = zeros(length(m.betas), q)
L[3, 1] = 1
F_test_sequence = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[3, "DenDF"], F_test_sequence[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[3, "FValue"], F_test_sequence[2], atol=1e-10, rtol=1e-5)

q = 3
L = zeros(length(m.betas), q)
L[4, 1] = 1
L[5, 2] = 1
L[6, 3] = 1
F_test_period = ftest_KR(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[1, "DenDF"], F_test_period[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[1, "FValue"], F_test_period[2], atol=1e-10, rtol=1e-4)

res = DataFrame(CSV.File("results/Results bioequivalence heterogeneous sas sw.csv"))

q = 1
L = zeros(length(m.betas), q)
L[2, 1] = 1
F_test_formulation = ftest_SW(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[2, "DenDF"], F_test_formulation[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[2, "FValue"], F_test_formulation[2], atol=1e-10, rtol=1e-7)

q = 1
L = zeros(length(m.betas), q)
L[3, 1] = 1
F_test_sequence = ftest_SW(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[3, "DenDF"], F_test_sequence[1], atol=1e-10, rtol=1e-6)
@test isapprox(res[3, "FValue"], F_test_sequence[2], atol=1e-10, rtol=1e-6)

q = 3
L = zeros(length(m.betas), q)
L[4, 1] = 1
L[5, 2] = 1
L[6, 3] = 1
F_test_period = ftest_SW(m, L; FIM_σ²=:observed_SAS_MATCHING)
@test isapprox(res[1, "DenDF"], F_test_period[1], atol=1e-10, rtol=1e-5)
@test isapprox(res[1, "FValue"], F_test_period[2], atol=1e-10, rtol=1e-7)
 =#
