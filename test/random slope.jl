using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results sleep study jmp.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(
    res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-6
)
@test isapprox(res[!, "DFDen"], kr.ν, atol=1e-10, rtol=1e-5)

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results sleep study sas kr.csv")))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-5, rtol=1e-10)
@test isapprox(res[!, "DF"], kr.ν, atol=1e-10, rtol=1e-5)

sw = small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results sleep study sas sw.csv")))
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-5, rtol=1e-10)
@test isapprox(res[!, "DF"], sw.ν, atol=1e-10, rtol=1e-5)

kr = small_sample_adjust(m, KenwardRoger(; fim=ExpectedFIM()))

res = DataFrame(CSV.File(joinpath(@__DIR__, "results", "Results sleep study lmertest.csv")))
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-10,
    rtol=1e-6,
)
@test isapprox(res[!, "coefficients.df"], kr.ν, atol=1e-10, rtol=1e-6)

sas_asycov_df = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study sas asycov.csv"))
)
W = vcov_varpar(m; fim=ObservedFIM())

row_residual_idx = findfirst(r -> r == "Residual", sas_asycov_df[!, "CovParm"])
row_int_idx = findfirst(
    i -> sas_asycov_df[i, "CovParm"] == "Intercept", 1:nrow(sas_asycov_df)
)
row_days_idx = findfirst(i -> sas_asycov_df[i, "CovParm"] == "days", 1:nrow(sas_asycov_df))

@test !isnothing(row_residual_idx)
@test !isnothing(row_int_idx)
@test !isnothing(row_days_idx)

sas_matrix = zeros(3, 3)
sas_matrix[1, 1] = sas_asycov_df[row_residual_idx, "CovP" * string(row_residual_idx)]
sas_matrix[1, 2] = sas_asycov_df[row_residual_idx, "CovP" * string(row_int_idx)]
sas_matrix[1, 3] = sas_asycov_df[row_residual_idx, "CovP" * string(row_days_idx)]

sas_matrix[2, 1] = sas_asycov_df[row_int_idx, "CovP" * string(row_residual_idx)]
sas_matrix[2, 2] = sas_asycov_df[row_int_idx, "CovP" * string(row_int_idx)]
sas_matrix[2, 3] = sas_asycov_df[row_int_idx, "CovP" * string(row_days_idx)]

sas_matrix[3, 1] = sas_asycov_df[row_days_idx, "CovP" * string(row_residual_idx)]
sas_matrix[3, 2] = sas_asycov_df[row_days_idx, "CovP" * string(row_int_idx)]
sas_matrix[3, 3] = sas_asycov_df[row_days_idx, "CovP" * string(row_days_idx)]

@test isapprox(W, sas_matrix, rtol=1e-4)

fm = @formula(reaction ~ 1 + days + (1 + days | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study corr sas kr.csv"))
)
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-5)
@test isapprox(res[!, "DF"], kr.ν, atol=1e-10, rtol=1e-4)

sas_asycov_df = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study corr sas asycov.csv"))
)
W = vcov_varpar(m; fim=ObservedFIM())
sas_names = sas_asycov_df[!, "CovParm"]
row_res = findfirst(==("Residual"), sas_names)
row_11 = findfirst(x -> occursin("UN(1,1)", x), sas_names)
row_21 = findfirst(x -> occursin("UN(2,1)", x), sas_names)
row_22 = findfirst(x -> occursin("UN(2,2)", x), sas_names)
@test !isnothing(row_res)
@test !isnothing(row_11)
@test !isnothing(row_21)
@test !isnothing(row_22)

inds = [row_res, row_11, row_21, row_22]
sas_matrix = zeros(4, 4)
for (i, r) in enumerate(inds)
    for (j, c) in enumerate(inds)
        val = sas_asycov_df[r, "CovP" * string(c)]
        sas_matrix[i, j] = val
    end
end
@test isapprox(W, sas_matrix, rtol=5e-4)

sw = small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study corr sas sw.csv"))
)
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-5, rtol=1e-5)
@test isapprox(res[!, "DF"], sw.ν, atol=1e-10, rtol=1e-4)

fm = @formula(reaction ~ 1 + days + days^2 + zerocorr(1 + days + days^2 | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study quadratic sas kr.csv"))
)
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-4)
@test isapprox(res[!, "DF"], kr.ν, atol=1e-10, rtol=0.5e-3)

sas_asycov_df = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study quadratic sas asycov.csv"))
)
W = vcov_varpar(m; fim=ObservedFIM())
row_res = findfirst(==("Residual"), sas_asycov_df[!, "CovParm"])
row_int = findfirst(x -> x == "Intercept", sas_asycov_df[!, "CovParm"])
row_d = findfirst(x -> x == "days", sas_asycov_df[!, "CovParm"])
row_d2 = findfirst(x -> occursin("days*days", x), sas_asycov_df[!, "CovParm"]) # SAS output often days*days for days^2
if isnothing(row_d2)
    row_d2 = findfirst(x -> occursin("days^2", x), sas_asycov_df[!, "CovParm"])
end
@test !isnothing(row_res)
@test !isnothing(row_int)
@test !isnothing(row_d)
@test !isnothing(row_d2)

inds = [row_res, row_int, row_d, row_d2]
sas_matrix = zeros(4, 4)
for (i, r) in enumerate(inds)
    for (j, c) in enumerate(inds)
        sas_matrix[i, j] = sas_asycov_df[r, "CovP" * string(c)]
    end
end
@test isapprox(W, sas_matrix, rtol=3e-4)

sw = small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study quadratic sas sw.csv"))
)
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-10, rtol=1e-4)
@test isapprox(res[!, "DF"], sw.ν, atol=1e-10, rtol=0.5e-3)

fm = @formula(reaction ~ 1 + days + days^2 + (1 + days + days^2 | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study corr quadratic sas kr.csv"))
)
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=0.5e-3)
@test isapprox(res[!, "DF"], kr.ν, atol=1e-10, rtol=1.5e-3)

sas_asycov_df = DataFrame(
    CSV.File(
        joinpath(@__DIR__, "results", "Results sleep study corr quadratic sas asycov.csv")
    ),
)
W = vcov_varpar(m; fim=ObservedFIM())
sas_names = sas_asycov_df[!, "CovParm"]
row_res = findfirst(==("Residual"), sas_names)
row_11 = findfirst(x -> occursin("UN(1,1)", x), sas_names)
row_21 = findfirst(x -> occursin("UN(2,1)", x), sas_names)
row_22 = findfirst(x -> occursin("UN(2,2)", x), sas_names)
row_31 = findfirst(x -> occursin("UN(3,1)", x), sas_names)
row_32 = findfirst(x -> occursin("UN(3,2)", x), sas_names)
row_33 = findfirst(x -> occursin("UN(3,3)", x), sas_names)

inds = [row_res, row_11, row_21, row_31, row_22, row_32, row_33]
@test !isnothing(row_res)
@test !isnothing(row_11)
@test !isnothing(row_21)
@test !isnothing(row_22)
@test !isnothing(row_31)
@test !isnothing(row_32)
@test !isnothing(row_33)

sas_matrix = zeros(7, 7)
for (i, r) in enumerate(inds)
    for (j, c) in enumerate(inds)
        sas_matrix[i, j] = sas_asycov_df[r, "CovP" * string(c)]
    end
end
@test isapprox(W, sas_matrix, rtol=3e-3)

sw = small_sample_adjust(m, Satterthwaite(; fim=ObservedFIM()))

res = DataFrame(
    CSV.File(joinpath(@__DIR__, "results", "Results sleep study corr quadratic sas sw.csv"))
)
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-10, rtol=0.5e-3)
@test isapprox(res[!, "DF"], sw.ν, atol=1e-10, rtol=1.5e-3)

fm = @formula(reaction ~ 1 + days + days^2 + (1 | subj) + (days + days^2 | subj))
m = fit(MixedModel, fm, df; REML=true)
kr = small_sample_adjust(m, KenwardRoger(; fim=ObservedFIM()))

res = DataFrame(
    CSV.File(
        joinpath(@__DIR__, "results", "Results sleep study some corr quadratic sas kr.csv")
    ),
)
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-10, rtol=1e-4)
@test isapprox(res[!, "DF"], kr.ν, atol=1e-10, rtol=1.0e-3)

sas_asycov_df = DataFrame(
    CSV.File(
        joinpath(
            @__DIR__, "results", "Results sleep study some corr quadratic sas asycov.csv"
        ),
    ),
)
W = vcov_varpar(m; fim=ObservedFIM())
sas_names = sas_asycov_df[!, "CovParm"]
row_res = findfirst(==("Residual"), sas_names)
row_int = findfirst(
    x -> x == "Intercept" || (occursin("Intercept", x) && !occursin("UN", x)), sas_names
)
row_un11 = findfirst(x -> occursin("UN(1,1)", x), sas_names)
row_un21 = findfirst(x -> occursin("UN(2,1)", x), sas_names)
row_un22 = findfirst(x -> occursin("UN(2,2)", x), sas_names)

inds = [5, 1, 2, 3, 4]

sas_matrix = zeros(5, 5)
for (i, r) in enumerate(inds)
    for (j, c) in enumerate(inds)
        sas_matrix[i, j] = sas_asycov_df[r, "CovP" * string(c)]
    end
end
@test isapprox(W, sas_matrix, rtol=2e-3)
