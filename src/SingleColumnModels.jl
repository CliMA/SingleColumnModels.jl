module SingleColumnModels

include(joinpath("Grids", "FiniteDifferenceGrids.jl"))
include(joinpath("StateVecs", "StateVecs.jl"))
include(joinpath("LinearSolvers", "TriDiagSolvers.jl"))
# include(joinpath("EDMF", "EDMF.jl"))

end # module
