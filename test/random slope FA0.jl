using CSV
using DataFrames
using Test

using MixedModelsSmallSample
using MixedModelsSmallSample.MixedModels
using MixedModelsSmallSample.LinearAlgebra: diag

# This test file validates FA0 (Cholesky) parameterization against SAS TYPE=FA0(2)
# For a 2-term random effect (intercept + days), FA0(2) is equivalent to full Cholesky

df = DataFrame(MixedModels.dataset(:sleepstudy))
fm = @formula(reaction ~ 1 + days + (1 + days | subj))
m = fit(MixedModel, fm, df; REML=true)

# Test that FA0 parameterization runs without error
@testset "FA0 parameterization executes" begin
    kr_fa0 = adjust_KR(m; FIM_σ²=:observed, parameterization=:FA0)
    sw_fa0 = adjust_SW(m; FIM_σ²=:observed, parameterization=:FA0)
    W_fa0 = MixedModelsSmallSample.vcov_varpar(m; FIM_σ²=:observed, parameterization=:FA0)

    @test size(W_fa0) == (4, 4)  # σ² + 3 Cholesky params (L11, L21, L22)
    @test length(kr_fa0.v) == 2
    @test length(sw_fa0.v) == 2
    @test all(isfinite.(kr_fa0.v))
    @test all(isfinite.(sw_fa0.v))
end

# Test that UN parameterization still works (backward compatibility)
@testset "UN parameterization backward compatibility" begin
    kr_un = adjust_KR(m; FIM_σ²=:observed, parameterization=:UN)
    sw_un = adjust_SW(m; FIM_σ²=:observed, parameterization=:UN)
    W_un = MixedModelsSmallSample.vcov_varpar(m; FIM_σ²=:observed, parameterization=:UN)

    @test size(W_un) == (4, 4)
    @test length(kr_un.v) == 2

    # Compare to SAS UN results (existing test data)
    res = DataFrame(
        CSV.File(joinpath(@__DIR__, "results", "Results sleep study corr sas kr.csv"))
    )
    @test isapprox(res[!, "Estimate"], kr_un.m.β, atol=1e-10, rtol=1e-10)
    @test isapprox(
        res[!, "StdErr"], sqrt.(diag(kr_un.varcovar_adjusted)), atol=1e-10, rtol=1e-5
    )
    @test isapprox(res[!, "DF"], kr_un.v, atol=1e-10, rtol=1e-4)
end

# Test FA0 against SAS FA0(2) results (when available)
# To be enabled once SAS results are provided
@testset "FA0 vs SAS FA0(2)" begin
    results_path = joinpath(@__DIR__, "results", "Results sleep study fa0 sas kr.csv")
    if isfile(results_path)
        res = DataFrame(CSV.File(results_path))
        kr_fa0 = adjust_KR(m; FIM_σ²=:observed, parameterization=:FA0)

        @test isapprox(res[!, "Estimate"], kr_fa0.m.β, atol=1e-10, rtol=1e-10)
        @test isapprox(
            res[!, "StdErr"], sqrt.(diag(kr_fa0.varcovar_adjusted)), atol=1e-5, rtol=1e-5
        )
        @test isapprox(res[!, "DF"], kr_fa0.v, atol=1e-2, rtol=1e-4)
    else
        @info "Skipping SAS FA0(2) comparison - results file not found: $results_path"
        @test_skip true
    end
end

# Test Satterthwaite FA0 (when SAS results available)
@testset "FA0 Satterthwaite vs SAS" begin
    results_path = joinpath(@__DIR__, "results", "Results sleep study fa0 sas sw.csv")
    if isfile(results_path)
        res = DataFrame(CSV.File(results_path))
        sw_fa0 = adjust_SW(m; FIM_σ²=:observed, parameterization=:FA0)

        @test isapprox(res[!, "Estimate"], sw_fa0.m.β, atol=1e-10, rtol=1e-10)
        @test isapprox(res[!, "StdErr"], sw_fa0.m.stderror, atol=1e-5, rtol=1e-5)
        @test isapprox(res[!, "DF"], sw_fa0.v, atol=1e-2, rtol=1e-4)
    else
        @info "Skipping SAS FA0(2) Satterthwaite comparison - results file not found: $results_path"
        @test_skip true
    end
end

# Test that default parameterization is :UN (backward compatibility)
@testset "Default parameterization is UN" begin
    kr_default = adjust_KR(m; FIM_σ²=:observed)
    kr_un = adjust_KR(m; FIM_σ²=:observed, parameterization=:UN)

    @test isapprox(kr_default.v, kr_un.v, rtol=1e-14)
    @test isapprox(kr_default.varcovar_adjusted, kr_un.varcovar_adjusted, rtol=1e-14)
end
