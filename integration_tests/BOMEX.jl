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
const test_data_dir = joinpath(pwd(), "output", "TestData")
const plot_dir = joinpath(scm_dir, "output", "bomex_edmf", "pycles_comparison")
mkpath(test_data_dir)

struct EarthParameterSet <: AbstractEarthParameterSet end
CLIMAParameters.Planet.MSLP(::EarthParameterSet) = 100000.0
const param_set = EarthParameterSet()

best_mse = OrderedDict()
best_mse[:q_tot_gm] = 6.8195421658167510e-01
best_mse[:a_up] = 1.2142815296731763e+02
best_mse[:w_up] = 7.9514345763404179e+01
best_mse[:q_tot_up] = 1.6612790419986577e+01
best_mse[:θ_liq_up] = 1.0041866755706111e+02
best_mse[:tke_en] = 9.6432737744200338e+01

# Define `compute_mse` and retrieve data file:
include(joinpath(@__DIR__, "utils", "compute_mse.jl"))
data_file = Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r")

@testset "BOMEX EDMF Solution Quality Assurance (QA) tests" begin

    grid, q, tmp, params = scm.EDMF.run(param_set, scm.EDMF.BOMEX())

    computed_mse = compute_mse(
        grid,
        q,
        tmp,
        params,
        data_file,
        "Bomex",
        best_mse,
        21600.0,
        plot_dir,
    )

    test_mse(computed_mse, best_mse, :q_tot_gm)
    test_mse(computed_mse, best_mse, :a_up)
    test_mse(computed_mse, best_mse, :w_up)
    test_mse(computed_mse, best_mse, :q_tot_up)
    test_mse(computed_mse, best_mse, :θ_liq_up)
    test_mse(computed_mse, best_mse, :tke_en)
end
