push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/SphericalTensors/src")
push!(LOAD_PATH, dirname(dirname(Base.@__DIR__)) * "/DMRG/src")
push!(LOAD_PATH, dirname(Base.@__DIR__) * "/src")

using Test, Random
using SphericalTensors, DMRG, Hamiltonians

Random.seed!(1342)

include("util.jl")

include("auxiliary.jl")

## algorithms
include("longrange.jl")


include("shortrange.jl")
