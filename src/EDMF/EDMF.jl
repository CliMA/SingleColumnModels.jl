"""
    Eddy-Diffusivity-Mass-Flux (EDMF)

Module for solving the sub-grid-scale equations
using the Eddy-Diffusivity-Mass-Flux (EDMF) model.
"""
module EDMF

using CLIMAParameters
using CLIMAParameters.Planet
using Thermodynamics
using UnPack
const TD = Thermodynamics
include("thermo_overload.jl")

using ..FiniteDifferenceGrids
using ..StateVecs
using ..TriDiagSolvers
using ..DomainDecompositions

include("cases.jl")
include("turb_conv_models.jl")
include("aux_funcs.jl")
include("init_params.jl")
include("updraft_vars.jl")
include("update_forcing.jl")
include("ref_states.jl")
include("initial_conditions.jl")
include("update_aux.jl")
include("apply_bcs.jl")
include("edmf_funcs.jl")
include("update_surface.jl")
include("update.jl")
include("main.jl")

end # module
