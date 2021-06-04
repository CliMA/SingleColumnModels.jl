using Test

# Unit tests / verification
for submodule in
    ["Grids", "DomainDecomposition", "StateVecs", "LinearSolvers", "PDEs"]

    println("Testing $submodule")
    include(joinpath(submodule, "runtests.jl"))
end

# Experiments / Integration tests
for submodule in ["BOMEX"]
    println("Testing $submodule")
    include(joinpath("..", "integration_tests", submodule * ".jl"))
end
