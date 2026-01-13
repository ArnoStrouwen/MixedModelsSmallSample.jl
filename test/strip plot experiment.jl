using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(CSV.File(joinpath(@__DIR__, "data", "Data battery cell Chapter 11.csv")))
rename!(df, "Whole Plots" => :WP)
rename!(df, "Subplots" => :SP)

fm = @formula(
    Y ~
        1 +
    (1 | WP) +
    (1 | SP) +
    X1 +
    X2 +
    X3 +
    X4 +
    X5 +
    X6 +
    (X2 + X3 + X4 + X5 + X6) & (X1) +
    (X3 + X4 + X5 + X6) & (X2) +
    (X4 + X5 + X6) & (X3) +
    (X5 + X6) & (X4) +
    (X6) & (X5)
)
m = fit(MixedModel, fm, df; REML=true)
kr = adjust_KR(m; FIM_σ²=:observed)

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results battery cell jmp.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-7, rtol=1e-7)
@test isapprox(
    res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-5, rtol=1e-10
)
@test isapprox(res[!, "DFDen"], kr.v, atol=1e-10, rtol=1e-5)

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results battery cell sas kr.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-6)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-4, rtol=1e-10)
@test isapprox(res[!, "DF"], kr.v, atol=1e-10, rtol=1e-4)

sw = adjust_SW(m; FIM_σ²=:observed)

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results battery cell sas sw.csv")))
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-6)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-4, rtol=1e-10)
@test isapprox(res[!, "DF"], sw.v, atol=1e-10, rtol=1e-4)

kr = adjust_KR(m; FIM_σ²=:expected)

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results battery cell lmertest.csv"))
)
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-10, rtol=1e-7)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-5,
    rtol=1e-10,
)
@test isapprox(res[!, "coefficients.df"], kr.v, atol=1e-10, rtol=1e-6)

sas_asycov_df = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results battery cell sas asycov.csv"))
)
W = MixedModelsSmallSample.vcov_varpar(m; FIM_σ²=:observed)

col_residual = findfirst(c -> c == "CovP1", names(sas_asycov_df))
col_sp = findfirst(c -> c == "CovP2", names(sas_asycov_df))
col_wp = findfirst(c -> c == "CovP3", names(sas_asycov_df))

row_residual_idx = findfirst(r -> r == "Residual", sas_asycov_df[!, "CovParm"])

row_wp_idx = findfirst(
    i -> sas_asycov_df[i, "CovParm"] == "Intercept" && sas_asycov_df[i, "Subject"] == "WP",
    1:nrow(sas_asycov_df),
)
row_sp_idx = findfirst(
    i -> sas_asycov_df[i, "CovParm"] == "Intercept" && sas_asycov_df[i, "Subject"] == "SP",
    1:nrow(sas_asycov_df),
)

@test !isnothing(col_residual)
@test !isnothing(col_sp)
@test !isnothing(col_wp)
@test !isnothing(row_residual_idx)
@test !isnothing(row_sp_idx)
@test !isnothing(row_wp_idx)

sas_matrix = zeros(3, 3)

sas_matrix[1, 1] = sas_asycov_df[row_residual_idx, "CovP" * string(row_residual_idx)]
sas_matrix[1, 2] = sas_asycov_df[row_residual_idx, "CovP" * string(row_wp_idx)]
sas_matrix[1, 3] = sas_asycov_df[row_residual_idx, "CovP" * string(row_sp_idx)]

sas_matrix[2, 1] = sas_asycov_df[row_wp_idx, "CovP" * string(row_residual_idx)]
sas_matrix[2, 2] = sas_asycov_df[row_wp_idx, "CovP" * string(row_wp_idx)]
sas_matrix[2, 3] = sas_asycov_df[row_wp_idx, "CovP" * string(row_sp_idx)]

sas_matrix[3, 1] = sas_asycov_df[row_sp_idx, "CovP" * string(row_residual_idx)]
sas_matrix[3, 2] = sas_asycov_df[row_sp_idx, "CovP" * string(row_wp_idx)]
sas_matrix[3, 3] = sas_asycov_df[row_sp_idx, "CovP" * string(row_sp_idx)]

@test isapprox(W, sas_matrix, rtol=1e-4)
