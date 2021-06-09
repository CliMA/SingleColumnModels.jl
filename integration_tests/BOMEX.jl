if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
using Test
using SingleColumnModels
using OrderedCollections
using SingleColumnModels.StateVecs
const scm = SingleColumnModels
using CLIMAParameters
using CLIMAParameters.Planet
import CLIMAParameters.Planet

ENV["GKSwstype"] = "nul"

const scm_dir = dirname(dirname(pathof(SingleColumnModels)));
output_dir = joinpath(scm_dir, "output", "BOMEX")
mkpath(output_dir)
const pycles_comparison = joinpath(output_dir, "pycles_comparison")
const plot_dir = joinpath(output_dir, "plots")

struct EarthParameterSet <: AbstractEarthParameterSet end
CLIMAParameters.Planet.MSLP(::EarthParameterSet) = 100000.0
const param_set = EarthParameterSet()

best_mse = OrderedDict()
best_mse[:q_tot_gm] = 6.8195420841336496e-01
best_mse[:a_up] = 1.2142815296731763e+02
best_mse[:w_up] = 7.9514302123362285e+01
best_mse[:q_tot_up] = 1.6612757475994467e+01
best_mse[:θ_liq_up] = 1.0041865638015472e+02
best_mse[:tke_en] = 9.6432752026034024e+01

# Define `compute_mse` and retrieve data file:
include(joinpath(@__DIR__, "utils", "compute_mse.jl"))
include(joinpath(@__DIR__, "utils", "plot_profiles.jl"))
data_file = Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r")

@testset "BOMEX EDMF Solution Quality Assurance (QA) tests" begin

    grid, q, tmp, params = scm.EDMF.run(param_set, scm.EDMF.BOMEX(), output_dir)
    computed_mse = compute_mse(
        grid,
        q,
        tmp,
        params,
        data_file,
        "Bomex",
        best_mse,
        21600.0,
        pycles_comparison,
    )

    # Export plots (must export before tests):
    skip_fields = ["a_gm", "u_en", "u_ud_1", "v_en", "v_gm", "v_ud_1", "w_gm"]
    prog = Dataset(joinpath(output_dir, "prog_vs_time.nc"), "r")
    plot_profiles(prog, plot_dir, skip_fields)
    # aux = Dataset(joinpath(output_dir, "aux_vs_time.nc"), "r")
    # plot_profiles(aux, plot_dir)

    test_mse(computed_mse, best_mse, :q_tot_gm)
    test_mse(computed_mse, best_mse, :a_up)
    test_mse(computed_mse, best_mse, :w_up)
    test_mse(computed_mse, best_mse, :q_tot_up)
    test_mse(computed_mse, best_mse, :θ_liq_up)
    test_mse(computed_mse, best_mse, :tke_en)
end
