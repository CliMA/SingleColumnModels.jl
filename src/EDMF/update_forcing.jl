#### Update forcing

####
#### Types
####

abstract type ForcingType end

struct NoForcing <: ForcingType end

struct StandardForcing{FT} <: ForcingType
  apply_subsidence::Bool
  apply_coriolis::Bool
  coriolis_param::FT
  StandardForcing(;apply_subsidence,apply_coriolis,coriolis_param::FT) where FT =
    new{FT}(apply_subsidence,apply_coriolis,coriolis_param)
end

####
#### Helper functions
####

function coriolis_force!(grid, q_tendencies, q, tmp, coriolis_param)
  gm, en, ud, sd, al = allcombinations(q)
  for k in over_elems_real(grid)
      q_tendencies[:u, k, gm] -= coriolis_param * (tmp[:vg, k] - q[:v, k, gm])
      q_tendencies[:v, k, gm] += coriolis_param * (tmp[:ug, k] - q[:u, k, gm])
  end
end

####
#### Update forcing methods
####

"""
    update_forcing!

Define forcing fields

 - `q_tendencies[:θ_liq, k, i]`
 - `q_tendencies[:q_tot, k, i]`

for all `k` and all `i`
"""
function update_forcing! end

function update_forcing!(q_tendencies::StateVec, tmp::StateVec, q::StateVec, grid::Grid, params, ::NoForcing) end

function update_forcing!(q_tendencies::StateVec, tmp::StateVec, q::StateVec, grid::Grid, params, forcing::StandardForcing)
  gm, en, ud, sd, al = allcombinations(q)
  @unpack params param_set

  for k in over_elems_real(grid)
    # Apply large-scale horizontal advection tendencies
    q_tendencies[:θ_liq, k, gm] += tmp[:dTdt, k]/exner_given_pressure(param_set, tmp[:p_0, k])
    q_tendencies[:q_tot, k, gm] += tmp[:dqtdt, k]
  end
  if forcing.apply_subsidence
    for k in over_elems_real(grid)
      # Apply large-scale subsidence tendencies
      q_tendencies[:θ_liq, k, gm] -= grad(q[:θ_liq, Dual(k), gm], grid) * tmp[:subsidence, k]
      q_tendencies[:q_tot, k, gm] -= grad(q[:q_tot, Dual(k), gm], grid) * tmp[:subsidence, k]
    end
  end

  if forcing.apply_coriolis
    coriolis_force!(grid, q_tendencies, q, tmp, forcing.coriolis_param)
  end
end
