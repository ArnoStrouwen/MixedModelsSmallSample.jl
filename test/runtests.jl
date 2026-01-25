using Test
using SafeTestsets

@testset "MixedModelsSmallSample.jl" begin
    @safetestset "Code quality (Aqua.jl)" begin
        using Aqua
        using MixedModelsSmallSample
        Aqua.test_all(MixedModelsSmallSample)
    end
    @safetestset "Code linting (JET.jl)" begin
        using JET
        using MixedModelsSmallSample
        JET.test_package(MixedModelsSmallSample; target_defined_modules=true)
    end
    @safetestset "blocked experiment" begin
        include("blocked experiment.jl")
    end
    @safetestset "split plot experiment" begin
        include("split plot experiment.jl")
    end
    @safetestset "strip plot experiment" begin
        include("strip plot experiment.jl")
    end
    @safetestset "random slope" begin
        include("random slope.jl")
    end
    @safetestset "categorical" begin
        include("bioequivalence homogeneous.jl")
    end
    @safetestset "categorical ML" begin
        include("bioequivalence homogeneous ml.jl")
    end
    @safetestset "heterogeneous bioequivalence" begin
        include("bioequivalence heterogeneous.jl")
    end
    @safetestset "varia" begin
        include("varia.jl")
    end
    @safetestset "random slope FA0" begin
        include("random slope FA0.jl")
    end
end
