using Test
using ASE
using JuLIP.Testing

h0("Testing ASE.jl")

@testset "ASE" begin
   @testset "Atoms"  begin include("testase.jl"); end
   @testset "Calculators" begin include("test_calculators.jl"); end
   @testset "JuLIP vs ASE" begin include("test_asevsjulip.jl"); end
end
