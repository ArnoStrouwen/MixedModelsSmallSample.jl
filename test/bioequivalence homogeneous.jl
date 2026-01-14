using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(CSV.File(joinpath(@__DIR__, "data", "Data bioequivalence.csv")))

fm = @formula(log_data ~ 1 + (1 | subject) + formulation + sequence + period)
# Match SAS reference levels: Period=4, Formulation=T, Sequence=BABA
m = fit(
    MixedModel,
    fm,
    df;
    REML=true,
    contrasts=Dict(
        :period => DummyCoding(; base=4),
        :formulation => DummyCoding(; base="T"),
        :sequence => DummyCoding(; base="BABA"),
    ),
)

kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res_pe = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results bioequivalence homogeneous sas kr.csv"))
)

julia_coefs = coef(m)
julia_se = stderror(m)
julia_dof = kr.ν

row_int = findfirst(r -> r["Effect"] == "Intercept", eachrow(res_pe))
@test isapprox(res_pe[row_int, "Estimate"], julia_coefs[1], atol=1e-5)
@test isapprox(res_pe[row_int, "StdErr"], julia_se[1], atol=1e-5)
@test isapprox(res_pe[row_int, "DF"], julia_dof[1], atol=1e-3)

row_form = findfirst(
    r -> r["Effect"] == "formulation" && r["formulation"] == "R", eachrow(res_pe)
)
@test isapprox(res_pe[row_form, "Estimate"], julia_coefs[2], atol=1e-5)
@test isapprox(res_pe[row_form, "StdErr"], julia_se[2], atol=1e-5)
@test isapprox(res_pe[row_form, "DF"], julia_dof[2], atol=1e-3)

row_seq = findfirst(
    r -> r["Effect"] == "sequence" && r["sequence"] == "ABAB", eachrow(res_pe)
)
@test isapprox(res_pe[row_seq, "Estimate"], julia_coefs[3], atol=1e-5)
@test isapprox(res_pe[row_seq, "StdErr"], julia_se[3], atol=1e-5)
@test isapprox(res_pe[row_seq, "DF"], julia_dof[3], atol=1e-3)

for (i, p) in enumerate([1, 2, 3])
    row_p = findfirst(
        r -> r["Effect"] == "period" && parse(Int, r["period"]) == p, eachrow(res_pe)
    )
    @test !isnothing(row_p)
    @test isapprox(res_pe[row_p, "Estimate"], julia_coefs[3 + i], atol=1e-5)
    @test isapprox(res_pe[row_p, "StdErr"], julia_se[3 + i], atol=1e-5)
    @test isapprox(res_pe[row_p, "DF"], julia_dof[3 + i], atol=1e-3)
end

res_tests = DataFrame(
    CSV.File(
        joinpath(
            @__DIR__, "results", "Results bioequivalence homogeneous sas tests3 kr.csv"
        ),
    ),
)

q = 1
L = zeros(length(m.betas), q)
L[2, 1] = 1
F_test_formulation = small_sample_ftest(m, L; method=KenwardRoger(; fim=ObservedFIM()))
row_form = findfirst(==("formulation"), res_tests[!, "Effect"])
@test isapprox(res_tests[row_form, "DenDF"], F_test_formulation[1], atol=1e-10, rtol=1e-5)
@test isapprox(res_tests[row_form, "FValue"], F_test_formulation[2], atol=1e-10, rtol=1e-4)

q = 1
L = zeros(length(m.betas), q)
L[3, 1] = 1
F_test_sequence = small_sample_ftest(m, L; method=KenwardRoger(; fim=ObservedFIM()))
row_seq = findfirst(==("sequence"), res_tests[!, "Effect"])
@test isapprox(res_tests[row_seq, "DenDF"], F_test_sequence[1], atol=1e-10, rtol=1e-5)
@test isapprox(res_tests[row_seq, "FValue"], F_test_sequence[2], atol=1e-10, rtol=1e-5)

q = 3
L = zeros(length(m.betas), q)
L[4, 1] = 1
L[5, 2] = 1
L[6, 3] = 1
F_test_period = small_sample_ftest(m, L; method=KenwardRoger(; fim=ObservedFIM()))
row_per = findfirst(==("period"), res_tests[!, "Effect"])
@test isapprox(res_tests[row_per, "NumDF"], q, atol=1e-10)
@test isapprox(res_tests[row_per, "DenDF"], F_test_period[1], atol=1e-10, rtol=1e-5)
@test isapprox(res_tests[row_per, "FValue"], F_test_period[2], atol=1e-10, rtol=1e-4)

sas_asycov_df = DataFrame(
    CSV.File(
        joinpath(@__DIR__, "results", "Results bioequivalence homogeneous sas asycov.csv")
    ),
)
W = vcov_varpar(m; fim=ObservedFIM())
row_res = findfirst(==("Residual"), sas_asycov_df[!, "CovParm"])
row_int = findfirst(x -> x == "Intercept", sas_asycov_df[!, "CovParm"])

inds = [row_res, row_int]
sas_matrix = zeros(2, 2)
for (i, r) in enumerate(inds)
    for (j, c) in enumerate(inds)
        sas_matrix[i, j] = sas_asycov_df[r, "CovP" * string(c)]
    end
end
@test isapprox(W, sas_matrix, rtol=1e-4)

sw = small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))

res_tests_sw = DataFrame(
    CSV.File(
        joinpath(
            @__DIR__, "results", "Results bioequivalence homogeneous sas tests3 sw.csv"
        ),
    ),
)

# Crucial test for bioequivalence
# Formulation (L vector construction)
# Julia Coefs: Int, Form:R, Seq:ABAB, Per:1, Per:2, Per:3
q = 1
L = zeros(length(m.betas), q)
L[2, 1] = 1
F_test_formulation = small_sample_ftest(m, L; method=Satterthwaite(; fim=ObservedFIM()))
row_form = findfirst(==("formulation"), res_tests_sw[!, "Effect"])
@test isapprox(
    res_tests_sw[row_form, "DenDF"], F_test_formulation[1], atol=1e-10, rtol=1e-5
)
@test isapprox(
    res_tests_sw[row_form, "FValue"], F_test_formulation[2], atol=1e-10, rtol=1e-6
)

q = 1
L = zeros(length(m.betas), q)
L[3, 1] = 1
F_test_sequence = small_sample_ftest(m, L; method=Satterthwaite(; fim=ObservedFIM()))
row_seq = findfirst(==("sequence"), res_tests_sw[!, "Effect"])
@test isapprox(res_tests_sw[row_seq, "DenDF"], F_test_sequence[1], atol=1e-10, rtol=1e-5)
@test isapprox(res_tests_sw[row_seq, "FValue"], F_test_sequence[2], atol=1e-10, rtol=1e-6)

q = 3
L = zeros(length(m.betas), q)
L[4, 1] = 1
L[5, 2] = 1
L[6, 3] = 1
F_test_period = small_sample_ftest(m, L; method=Satterthwaite(; fim=ObservedFIM()))
row_per = findfirst(==("period"), res_tests_sw[!, "Effect"])
@test isapprox(res_tests_sw[row_per, "DenDF"], F_test_period[1], atol=1e-10, rtol=1e-5)
@test isapprox(res_tests_sw[row_per, "FValue"], F_test_period[2], atol=1e-10, rtol=1e-6)
