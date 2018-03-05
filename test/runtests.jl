
using ASE
using Base.Test

println("-------------------")
println(" Testing ASE.jl")
println("-------------------")

@testset "ASE"  begin
include("testase.jl")
end
