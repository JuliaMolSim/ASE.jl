using Test
using ASE
include("aux.jl")

h0("Testing ASE.jl")

@testset "ASE" begin
   # @testset "Build" begin include("test_build.jl"); end
   @testset "Atoms"  begin include("testase.jl"); end
   # @testset "Calculators" begin include("test_calculators.jl"); end
end
