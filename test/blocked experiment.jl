using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

df = DataFrame(CSV.File("data/Data Pastry Dough Experiment Chapter 7.csv"))
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
kr = adjust_KR(m; FIM_σ²=:observed_SAS_MATCHING)

res = DataFrame(CSV.File("results/Results pastry dough jmp.csv"))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-9, rtol=1e-9)
@test isapprox(
    res[!, "Std Error"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-9, rtol=1e-10
)
@test isapprox(res[!, "DFDen"], kr.v, atol=1e-10, rtol=1e-8)

res = DataFrame(CSV.File("results/Results pastry dough sas kr.csv"))
@test isapprox(res[!, "Estimate"], kr.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sqrt.(diag(kr.varcovar_adjusted)), atol=1e-9, rtol=1e-10)
@test isapprox(res[!, "DF"], kr.v, atol=1e-10, rtol=1e-8)

sw = adjust_SW(m; FIM_σ²=:observed_SAS_MATCHING)

res = DataFrame(CSV.File("results/Results pastry dough sas sw.csv"))
@test isapprox(res[!, "Estimate"], sw.m.β, atol=1e-10, rtol=1e-10)
@test isapprox(res[!, "StdErr"], sw.m.stderror, atol=1e-9, rtol=1e-10)
@test isapprox(res[!, "DF"], sw.v, atol=1e-10, rtol=1e-8)

kr = adjust_KR(m; FIM_σ²=:expected)

res = DataFrame(CSV.File("results/Results pastry dough lmertest.csv"))
res = vcat(res, res[5:7, :])
deleteat!(res, 5:7)
@test isapprox(res[!, "coefficients.Estimate"], kr.m.β, atol=1e-8, rtol=1e-8)
@test isapprox(
    res[!, "coefficients.Std..Error"],
    sqrt.(diag(kr.varcovar_adjusted)),
    atol=1e-6,
    rtol=1e-7,
)
@test isapprox(res[!, "coefficients.df"], kr.v, atol=1e-7, rtol=1e-7)
