###### Pressure models

abstract type PressureModel end

struct SCAMPyPressure{FT} <: PressureModel
  buoy_coeff::FT
  drag_coeff::FT
  plume_spacing::FT
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

