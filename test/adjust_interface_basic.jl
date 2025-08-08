using Test
using MixedModelsSmallSample

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

# Note: Full functional tests requiring MixedModels data are included 
# in other test files that have proper test environments set up.
