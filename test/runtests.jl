
using ASE
using Base.Test

println("-------------------")
println(" Testing ASE.jl")
println("-------------------")

@testset "ASE" begin
   @testset "Build" begin include("test_build.jl"); end
   @testset "Atoms"  begin include("testase.jl"); end
   @testset "Calculators" begin include("test_calculators.jl"); end
end
