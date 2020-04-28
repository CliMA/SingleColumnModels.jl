using Test

using SingleColumnModels
using SingleColumnModels.StateVecs

using CLIMAParameters
import CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
CLIMAParameters.Planet.MSLP(::EarthParameterSet) = 100000.0

using CLIMAParameters.Planet
const param_set = EarthParameterSet()

const test_data_dir = joinpath(pwd(),"output","TestData")
mkpath(test_data_dir)

#### Accepting a new solution
# After a careful review of ALL of the solution results,
# set `accept_new_solution = true`, and run once. Then,
# reset to `false` before committing. This parameter should
# never be `true` for commits.
const accept_new_solution = true # - true always pass tests

scms = SingleColumnModels
@testset "Integration test: EDMF equations (BOMEX)" begin
  @show param_set
  grid, q, tmp = scms.EDMF.run(param_set, scms.EDMF.BOMEX())

  if accept_new_solution || !isfile(joinpath(test_data_dir, "q_expected.csv"))
    export_state(q, grid, test_data_dir, "q_expected.csv")
  end
  if accept_new_solution || !isfile(joinpath(test_data_dir, "tmp_expected.csv"))
    export_state(tmp, grid, test_data_dir, "tmp_expected.csv")
  end

  gm, en, ud, sd, al = allcombinations(q)
  FT = eltype(grid)
  q_expected = deepcopy(q)
  tmp_expected = deepcopy(tmp)
  assign!(q_expected, grid, FT(0))
  assign!(tmp_expected, grid, FT(0))

  import_state!(q_expected, grid, test_data_dir, "q_expected.csv")
  import_state!(tmp_expected, grid, test_data_dir, "tmp_expected.csv")

  D_q = compare(q, q_expected, grid, eps(Float32))
  D_tmp = compare(tmp, tmp_expected, grid, eps(Float32))
  # plot_state(q, grid, test_data_dir, "q")

  @test all(D_q[:a])
  @test all(D_q[:w])
  @test all(D_q[:θ_liq])
  @test all(D_q[:q_tot])
  @test all(D_q[:tke])
  @test all(D_tmp[:q_liq])

end

