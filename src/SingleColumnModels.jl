module SingleColumnModels

include(joinpath("Grids", "FiniteDifferenceGrids.jl"))
include(joinpath("DomainDecompositions", "DomainDecompositions.jl"))
include(joinpath("StateVecs", "StateVecs.jl"))
include(joinpath("EDMF", "EDMF.jl"))

end # module
