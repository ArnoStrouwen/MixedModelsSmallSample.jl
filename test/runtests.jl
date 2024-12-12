using KenwardRoger
using Test
using Aqua
using JET

@testset "KenwardRoger.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(KenwardRoger)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(KenwardRoger; target_defined_modules = true)
    end
    @testset "Scientific tests" begin
        include("kenward_roger.jl")
    end
end
