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
output_dir = joinpath(scm_dir, "output", "DYCOMS")
mkpath(output_dir)
const pycles_comparison = joinpath(output_dir, "pycles_comparison")
const plot_dir = joinpath(output_dir, "plots")

struct EarthParameterSet <: AbstractEarthParameterSet end
CLIMAParameters.Planet.MSLP(::EarthParameterSet) = 100000.0
CLIMAParameters.Planet.kappa_d(::EarthParameterSet) = 0.2827586206896552
const param_set = EarthParameterSet()

best_mse = OrderedDict()
best_mse[:q_tot_gm] = 8.6933572671880666e-01
best_mse[:a_up] = 1.3643220197659051e+02
best_mse[:w_up] = 4.9458009107657077e+01
best_mse[:q_tot_up] = 1.3778729479726136e+01
best_mse[:θ_liq_up] = 1.0043173075802038e+02
best_mse[:tke_en] = 9.2145436040979149e+01

# Define `compute_mse` and retrieve data file:
include(joinpath(scm_dir, "integration_tests", "utils", "compute_mse.jl"))
include(joinpath(scm_dir, "integration_tests", "utils", "plot_profiles.jl"))
data_file = Dataset(joinpath(PyCLES_output_dataset_path, "DYCOMS_RF01.nc"), "r")

@testset "DYCOMS EDMF Solution Quality Assurance (QA) tests" begin

    grid, q, aux, params =
        scm.EDMF.run(param_set, scm.EDMF.DYCOMS(), output_dir)
    computed_mse = compute_mse(
        grid,
        q,
        aux,
        params,
        data_file,
        "Dycoms",
        best_mse,
        21600.0,
        pycles_comparison,
    )

    # Export plots (must export before tests):
    skip_fields = ["a_gm", "u_en", "u_ud_1", "v_en", "v_gm", "v_ud_1", "w_gm"]
    prog = Dataset(joinpath(output_dir, "prog_vs_time.nc"), "r")
    plot_profiles(prog, plot_dir; skip_fields = skip_fields)
    # aux = Dataset(joinpath(output_dir, "aux_vs_time.nc"), "r")
    # plot_profiles(aux, plot_dir)

    test_mse(computed_mse, best_mse, :q_tot_gm)
    test_mse(computed_mse, best_mse, :a_up)
    test_mse(computed_mse, best_mse, :w_up)
    test_mse(computed_mse, best_mse, :q_tot_up)
    test_mse(computed_mse, best_mse, :θ_liq_up)
    test_mse(computed_mse, best_mse, :tke_en)
end
