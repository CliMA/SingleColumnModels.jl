#### UpdateForcing

function update_forcing! end

function update_forcing!(q_tendencies::StateVec, tmp::StateVec, q::StateVec, grid::Grid, params, ::BOMEX)
  gm, en, ud, sd, al = allcombinations(DomainIdx(q))
  @unpack params param_set
  z = grid.zc
  for k in over_elems_real(grid)
    # Geostrophic velocity profiles. vg = 0
    tmp[:ug, k, gm] = -10.0 + (1.8e-3)*z[k]
    # Set large-scale cooling
    if z[k] <= 1500.0
      tmp[:dTdt, k, gm] =  (-2.0/(3600 * 24.0))  * exner_given_pressure(param_set, tmp[:p_0, k])
    else
      tmp[:dTdt, k, gm] = (-2.0/(3600 * 24.0) + (z[k] - 1500.0)
                          * (0.0 - -2.0/(3600 * 24.0)) / (3000.0 - 1500.0)) * exner_given_pressure(param_set, tmp[:p_0, k])
    end
    # Set large-scale drying
    if z[k] <= 300.0
      tmp[:dqtdt, k, gm] = -1.2e-8   #kg/(kg * s)
    end
    if z[k] > 300.0 and z[k] <= 500.0
      tmp[:dqtdt, k, gm] = -1.2e-8 + (z[k] - 300.0)*(0.0 - -1.2e-8)/(500.0 - 300.0) #kg/(kg * s)
    end

    #Set large scale subsidence
    if z[k] <= 1500.0
      tmp[:subsidence, k, gm] = 0.0 + z[k]*(-0.65/100.0 - 0.0)/(1500.0 - 0.0)
    end
    if z[k] > 1500.0 and z[k] <= 2100.0
      tmp[:subsidence, k, gm] = -0.65/100 + (z[k] - 1500.0)* (0.0 - -0.65/100.0)/(2100.0 - 1500.0)
    end
  end
end


function coriolis_force!(grid, q_tendencies, q, tmp, coriolis_param)
  gm, en, ud, sd, al = allcombinations(q)
  @inbounds for k in over_elems_real(grid)
    q_tendencies[:u, k, gm] -= coriolis_param * (tmp[:vg, k] - q[:v, k, gm])
    q_tendencies[:v, k, gm] += coriolis_param * (tmp[:ug, k] - q[:u, k, gm])
  end
end

abstract type ForcingType end
struct NoForcing <: ForcingType end
struct StandardForcing{FT} <: ForcingType
  apply_subsidence::Bool
  apply_coriolis::Bool
  coriolis_param::FT
  StandardForcing(;apply_subsidence,apply_coriolis,coriolis_param::FT) where FT =
    new{FT}(apply_subsidence,apply_coriolis,coriolis_param)
end

function update_forcing!(q_tendencies::StateVec, tmp::StateVec, q::StateVec, grid::Grid, params, ::ForcingType)
end

function update_forcing!(q_tendencies::StateVec, tmp::StateVec, q::StateVec, grid::Grid, params, model::StandardForcing)
  gm, en, ud, sd, al = allcombinations(q)
  @unpack params param_set

  for k in over_elems_real(grid)
    # Apply large-scale horizontal advection tendencies
    q_tendencies[:θ_liq, k, gm] += tmp[:dTdt, k]/exner_given_pressure(param_set, tmp[:p_0, k])
    q_tendencies[:q_tot, k, gm] += tmp[:dqtdt, k]
  end
  if model.apply_subsidence
    for k in over_elems_real(grid)
      # Apply large-scale subsidence tendencies
      q_tendencies[:θ_liq, k, gm] -= grad(q[:θ_liq, Dual(k), gm], grid) * tmp[:subsidence, k]
      q_tendencies[:q_tot, k, gm] -= grad(q[:q_tot, Dual(k), gm], grid) * tmp[:subsidence, k]
    end
  end

  if model.apply_coriolis
    coriolis_force!(grid, q_tendencies, q, tmp, model.coriolis_param)
  end

end
