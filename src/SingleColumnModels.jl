module SingleColumnModels

import ClimateMachineCore.DataLayouts: VF
using ClimateMachineCore: Fields, Spaces

export FaceField, CentField

function CentField(space::Spaces.FaceFiniteDifferenceSpace, vals::NamedTuple)
    S = typeof(vals)
    FT = Spaces.undertype(space)
    n_fields = Int(sizeof(vals)/sizeof(FT))
    array = zeros(FT, Spaces.n_cells(space), n_fields)
    return Fields.Field(VF{S}(array), space)
end

function FaceField(space::Spaces.FaceFiniteDifferenceSpace, vals::NamedTuple)
    S = typeof(vals)
    FT = Spaces.undertype(space)
    n_fields = Int(sizeof(vals)/sizeof(FT))
    array = zeros(FT, Spaces.n_faces(space), n_fields)
    return Fields.Field(VF{S}(array), space)
end


include(joinpath("Grids", "FiniteDifferenceGrids.jl"))
include(joinpath("DomainDecompositions", "DomainDecompositions.jl"))
include(joinpath("StateVecs", "StateVecs.jl"))
include(joinpath("LinearSolvers", "TriDiagSolvers.jl"))
include(joinpath("EDMF", "EDMF.jl"))

end # module
