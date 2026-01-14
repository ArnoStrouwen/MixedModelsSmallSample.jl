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
    kr_fa0 = small_sample_adjust(
        m, KenwardRoger(fim=ObservedFIM(), parameterization=FactorAnalytic())
    )
    sw_fa0 = small_sample_adjust(
        m, Satterthwaite(fim=ObservedFIM(), parameterization=FactorAnalytic())
    )
    W_fa0 = vcov_varpar(m; fim=ObservedFIM(), parameterization=FactorAnalytic())

    @test size(W_fa0) == (4, 4)  # σ² + 3 Cholesky params (L11, L21, L22)
    @test length(kr_fa0.ν) == 2
    @test length(sw_fa0.ν) == 2
    @test all(isfinite.(kr_fa0.ν))
    @test all(isfinite.(sw_fa0.ν))
end

# Test that UN parameterization still works (backward compatibility)
@testset "UN parameterization backward compatibility" begin
    kr_un = small_sample_adjust(
        m, KenwardRoger(fim=ObservedFIM(), parameterization=Unstructured())
    )
    sw_un = small_sample_adjust(
        m, Satterthwaite(fim=ObservedFIM(), parameterization=Unstructured())
    )
    W_un = vcov_varpar(m; fim=ObservedFIM(), parameterization=Unstructured())

    @test size(W_un) == (4, 4)
    @test length(kr_un.ν) == 2

    # Compare to SAS UN results (existing test data)
    res = DataFrame(
        CSV.File(joinpath(@__DIR__, "results", "Results sleep study corr sas kr.csv"))
    )
    @test isapprox(res[!, "Estimate"], kr_un.m.β, atol=1e-10, rtol=1e-10)
    @test isapprox(
        res[!, "StdErr"], sqrt.(diag(kr_un.varcovar_adjusted)), atol=1e-10, rtol=1e-5
    )
    @test isapprox(res[!, "DF"], kr_un.ν, atol=1e-10, rtol=1e-4)
end

# Test FA0 against SAS FA0(2) results (when available)
# To be enabled once SAS results are provided
@testset "FA0 vs SAS FA0(2)" begin
    results_path = joinpath(@__DIR__, "results", "Results sleep study fa0 sas kr.csv")
    if isfile(results_path)
        res = DataFrame(CSV.File(results_path))
        kr_fa0 = small_sample_adjust(
            m, KenwardRoger(fim=ObservedFIM(), parameterization=FactorAnalytic())
        )

        @test isapprox(res[!, "Estimate"], kr_fa0.m.β, atol=1e-10, rtol=1e-10)
        @test isapprox(
            res[!, "StdErr"], sqrt.(diag(kr_fa0.varcovar_adjusted)), atol=1e-5, rtol=1e-5
        )
        @test isapprox(res[!, "DF"], kr_fa0.ν, atol=1e-2, rtol=1e-4)
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
        sw_fa0 = small_sample_adjust(
            m, Satterthwaite(fim=ObservedFIM(), parameterization=FactorAnalytic())
        )

        @test isapprox(res[!, "Estimate"], sw_fa0.m.β, atol=1e-10, rtol=1e-10)
        @test isapprox(res[!, "StdErr"], sw_fa0.m.stderror, atol=1e-5, rtol=1e-5)
        @test isapprox(res[!, "DF"], sw_fa0.ν, atol=1e-2, rtol=1e-4)
    else
        @info "Skipping SAS FA0(2) Satterthwaite comparison - results file not found: $results_path"
        @test_skip true
    end
end

# Test that default parameterization is :UN (backward compatibility)
@testset "Default parameterization is UN" begin
    kr_default = small_sample_adjust(m, KenwardRoger(fim=ObservedFIM()))
    kr_un = small_sample_adjust(
        m, KenwardRoger(fim=ObservedFIM(), parameterization=Unstructured())
    )

    @test isapprox(kr_default.ν, kr_un.ν, rtol=1e-14)
    @test isapprox(kr_default.varcovar_adjusted, kr_un.varcovar_adjusted, rtol=1e-14)
end
