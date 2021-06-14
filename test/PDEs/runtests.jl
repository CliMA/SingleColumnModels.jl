using Test

using SingleColumnModels.FiniteDifferenceGrids
using SingleColumnModels.StateVecs
using LinearAlgebra

output_root = joinpath("..", "output", "tests", "PDEs")

include("HeatEquation.jl")
include("AdvectionEquation.jl")
