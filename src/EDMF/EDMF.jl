"""
    Eddy-Diffusivity-Mass-Flux (EDMF)

Module for solving the sub-grid-scale equations
using the Eddy-Diffusivity-Mass-Flux (EDMF) model.
"""
module EDMF

using CLIMAParameters
using CLIMAParameters.Planet
using MoistThermodynamics
include("moist_thermo_overload.jl")

using ..FiniteDifferenceGrids
using ..StateVecs
using ..TriDiagSolvers

include("utilities.jl")
include("cases.jl")
include("dir_trees.jl")
include("turb_conv_models.jl")
include("aux_funcs.jl")
include("init_params.jl")
include("updraft_vars.jl")
include("forcing_funcs.jl")
include("ref_states.jl")
include("initial_conditions.jl")
include("update_aux.jl")
include("apply_bcs.jl")
include("edmf_funcs.jl")
include("surface_funcs.jl")
include("update.jl")
include("process_results.jl")
include("main.jl")

end # module
