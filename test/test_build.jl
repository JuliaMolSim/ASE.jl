using PyCall

include("../deps/build.jl")

println("Test Dependency Build")
@test pip("pip") == nothing