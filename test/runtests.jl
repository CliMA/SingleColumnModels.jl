using Test

# Unit tests / verification
for submodule in ["Grids", "DomainDecomposition", "StateVecs", "PDEs"]
    println("Testing $submodule")
    include(joinpath(submodule, "runtests.jl"))
end
