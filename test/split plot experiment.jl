using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(CSV.File(joinpath(@__DIR__, "data", "Data wind tunnel Chapter 10.csv")))
rename!(df, "Whole Plots" => :WP)

fm = @formula(
    EFFICIENCY ~
        1 +
    (1 | WP) +
    FRH +
    RRH +
    YA +
    GC +
    FRH & RRH +
    FRH & YA +
    FRH & GC +
    RRH & YA +
    RRH & GC +
    YA & GC +
    FRH & FRH +
    RRH & RRH +
    YA & YA +
    GC & GC
)
m = fit(MixedModel, fm, df; REML=true)

kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results wind tunnel jmp.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-9, rtol=1e-8)
@test isapprox(
    res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-5
)
@test isapprox(res[!, "DFDen"], kr.ν, atol=1e-10, rtol=1e-6)

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results wind tunnel sas kr.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-9, rtol=1e-8)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-7, rtol=1e-10)
@test isapprox(res[!, "DF"], kr.ν, atol=1e-10, rtol=1e-6)

sw = small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results wind tunnel sas sw.csv")))
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-9, rtol=1e-8)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-7, rtol=1e-10)
@test isapprox(res[!, "DF"], sw.ν, atol=1e-10, rtol=1e-6)

kr = small_sample_adjust(m, KenwardRoger(; fim=ExpectedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results wind tunnel lmertest.csv")))
res = vcat(res, res[6:9, :])
deleteat!(res, 6:9)
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-10, rtol=1e-8)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-6,
    rtol=1e-9,
)
@test isapprox(res[!, "coefficients.df"], kr.ν, atol=1e-9, rtol=1e-7)

sas_asycov_df = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results wind tunnel sas asycov.csv"))
)
W = vcov_varpar(m; fim=ObservedFIM())

col_residual = findfirst(c -> c == "CovP1", names(sas_asycov_df))
col_intercept = findfirst(c -> c == "CovP2", names(sas_asycov_df))

row_residual_idx = findfirst(r -> occursin("Residual", r), sas_asycov_df[!, "CovParm"])
row_intercept_idx = findfirst(r -> occursin("Intercept", r), sas_asycov_df[!, "CovParm"])

@test !isnothing(col_residual)
@test !isnothing(col_intercept)
@test !isnothing(row_residual_idx)
@test !isnothing(row_intercept_idx)

sas_matrix = zeros(2, 2)

sas_matrix[1, 1] = sas_asycov_df[row_residual_idx, "CovP" * string(row_residual_idx)]
sas_matrix[1, 2] = sas_asycov_df[row_residual_idx, "CovP" * string(row_intercept_idx)]
sas_matrix[2, 1] = sas_asycov_df[row_intercept_idx, "CovP" * string(row_residual_idx)]
sas_matrix[2, 2] = sas_asycov_df[row_intercept_idx, "CovP" * string(row_intercept_idx)]

@test isapprox(W, sas_matrix, rtol=1e-4)
