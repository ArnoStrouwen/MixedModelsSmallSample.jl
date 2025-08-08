using Test
using MixedModelsSmallSample

# Import dependencies for functional tests
using MixedModels
using DataFrames

# Test that the types are properly defined
@testset "Abstract type hierarchy" begin
    @test KenwardRoger <: AbstractLMMSS
    @test Satterthwaite <: AbstractLMMSS

    # Test construction
    kr = KenwardRoger()
    sw = Satterthwaite()
    @test kr isa KenwardRoger
    @test sw isa Satterthwaite
    @test kr isa AbstractLMMSS
    @test sw isa AbstractLMMSS
end

# Test the new adjust function with sample data if available
@testset "adjust function interface" begin
    # Try to test with real data, but gracefully handle network failures
    try
        df = DataFrame(MixedModels.dataset(:sleepstudy))
        fm = @formula(reaction ~ 1 + days + zerocorr(1 + days | subj))
        m = fit(MixedModel, fm, df; REML=true)  # Use REML=true

        # Test default method (should be KenwardRoger)
        result_default = adjust(m)
        @test result_default isa LinearMixedModelKR

        # Test explicit KenwardRoger method  
        result_kr = adjust(m; method=KenwardRoger())
        @test result_kr isa LinearMixedModelKR
        @test result_kr.m === m  # should contain original model

        # Test Satterthwaite method
        result_sw = adjust(m; method=Satterthwaite())
        @test result_sw isa LinearMixedModelSW
        @test result_sw.m === m  # should contain original model

        # Test that results are equivalent to original functions
        kr_old = adjust_KR(m)
        sw_old = adjust_SW(m)
        @test result_kr.varcovar_adjusted ≈ kr_old.varcovar_adjusted
        @test result_sw.v ≈ sw_old.v

        println("Full functional tests passed!")

    catch e
        if isa(e, SystemError) ||
            contains(string(e), "network") ||
            contains(string(e), "download")
            @warn "Skipping functional tests due to network restrictions: $e"
        else
            rethrow(e)  # Re-throw unexpected errors
        end
    end
end
