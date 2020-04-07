###### Pressure models

abstract type PressureModel end

Base.@kwdef struct SCAMPyPressure{FT} <: PressureModel
  buoy_coeff::FT=FT(0)
  drag_coeff::FT=FT(0)
  plume_spacing::FT=FT(0)
end

function compute_pressure!(grid::Grid{FT}, q, tmp, params, model::SCAMPyPressure) where FT
  gm, en, ud, sd, al = allcombinations(q)
  p_coeff = model.drag_coeff/model.plume_spacing

  @inbounds for i in ud
    @inbounds for k in over_elems_real(grid)
      a_k = q[:a, k, i]
      w_env = q[:w, k, en]
      ρ_k = tmp[:ρ_0, k]
      w_i = q[:w, k, i]
      B_k = tmp[:buoy, k, i]
      ρa_k = ρ_k * a_k

      press_buoy = - ρa_k * B_k * model.buoy_coeff
      press_drag = - ρa_k * (p_coeff * (w_i - w_env)^2/sqrt(a_k))
      tmp[:nh_press, k, i] = press_buoy + press_drag
    end
  end
end

function compute_tke_pressure!(grid::Grid{FT}, q, tmp, tmp_O2, cv, params, model::SCAMPyPressure) where FT
  gm, en, ud, sd, al = allcombinations(q)
  p_coeff = model.drag_coeff/model.plume_spacing
  @inbounds for k in over_elems_real(grid)
    tmp_O2[cv][:press, k] = FT(0)
    @inbounds for i in ud
      wu_half = q[:w, k, i]
      we_half = q[:w, k, en]
      a_i = q[:a, k, i]
      ρ_0_k = tmp[:ρ_0, k]
      press_buoy = -1 * ρ_0_k * a_i * tmp[:buoy, k, i] * model.buoy_coeff
      press_drag_coeff = -1 * ρ_0_k * sqrt(a_i) * p_coeff
      press_drag = press_drag_coeff * (wu_half - we_half)*abs(wu_half - we_half)
      tmp_O2[cv][:press, k] += (we_half - wu_half) * (press_buoy + press_drag)
    end
  end
end
