using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(
    CSV.File(joinpath(@__DIR__, "data", "Data Pastry Dough Experiment Chapter 7.csv"))
)
rename!(df, "Flow Rate" => :FR, "Moisture Content" => :MC, "Screw Speed" => :SS)
rename!(df, "Longitudinal Expansion Index" => :LEI)
function one_minus_one_coding!(x)
    minx = minimum(x)
    maxx = maximum(x)
    span = maxx - minx
    x .-= minx .+ span / 2
    x ./= (span / 2)
    return nothing
end
one_minus_one_coding!(df[!, :FR])
one_minus_one_coding!(df[!, :MC])
one_minus_one_coding!(df[!, :SS])
fm = @formula(
    LEI ~
        1 +
    (1 | Day) +
    FR +
    MC +
    SS +
    FR & MC +
    FR & SS +
    MC & SS +
    FR & FR +
    MC & MC +
    SS & SS
)
m = fit(MixedModel, fm, df; REML=true)
kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results pastry dough jmp.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-9, rtol=1e-9)
@test isapprox(res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-9, rtol=1e-8)
@test isapprox(res[!, "DFDen"], kr.ν, atol=1e-10, rtol=1e-7)

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results pastry dough sas kr.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-9, rtol=1e-8)
@test isapprox(res[!, "DF"], kr.ν, atol=1e-10, rtol=1e-7)

sw = small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results pastry dough sas sw.csv")))
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-9, rtol=1e-8)
@test isapprox(res[!, "DF"], sw.ν, atol=1e-10, rtol=1e-7)

kr = small_sample_adjust(m, KenwardRoger(; fim=ExpectedFIM()))

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results pastry dough lmertest.csv"))
)
res = vcat(res, res[5:7, :])
deleteat!(res, 5:7)
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-8, rtol=1e-8)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-6,
    rtol=1e-7,
)
@test isapprox(res[!, "coefficients.df"], kr.ν, atol=1e-7, rtol=1e-7)

sas_asycov_df = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results pastry dough sas asycov.csv"))
)
W = vcov_varpar(m; fim=ObservedFIM())
col_intercept = findfirst(c -> c == "CovP2", names(sas_asycov_df))
col_residual = findfirst(c -> c == "CovP1", names(sas_asycov_df))

row_residual_idx = findfirst(r -> occursin("Residual", r), sas_asycov_df[!, "CovParm"])
row_intercept_idx = findfirst(r -> occursin("Intercept", r), sas_asycov_df[!, "CovParm"])

sas_matrix = zeros(2, 2)

sas_matrix[1, 1] = sas_asycov_df[row_residual_idx, "CovP" * string(row_residual_idx)]
sas_matrix[1, 2] = sas_asycov_df[row_residual_idx, "CovP" * string(row_intercept_idx)]
sas_matrix[2, 1] = sas_asycov_df[row_intercept_idx, "CovP" * string(row_residual_idx)]
sas_matrix[2, 2] = sas_asycov_df[row_intercept_idx, "CovP" * string(row_intercept_idx)]

@test isapprox(W, sas_matrix, rtol=1e-7)
