###### Pressure models

abstract type PressureModel end

Base.@kwdef struct SCAMPyPressure{FT} <: PressureModel
  buoy_coeff::FT=FT(0)
  drag_coeff::FT=FT(0)
  plume_spacing::FT=FT(0)
end

Base.@kwdef struct NormalMode{FT} <: PressureModel
  nh_avd_coef::FT=FT(0)
  buoy_coef::FT=FT(0)
  nh_drag_coef::FT=FT(0)
end

"""
    compute_pressure!

Define pressure field

 - `tmp[:nh_press, k, i]`

for all `k` and all `i`
"""
function compute_pressure! end

function compute_pressure!(grid::Grid{FT}, UpdVar, q, tmp, params, model::SCAMPyPressure) where FT
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

function compute_pressure!(grid::Grid{FT}, UpdVar, q, tmp, params, model::NormalMode) where FT
  gm, en, ud, sd, al = allcombinations(q)

  @inbounds for i in ud
    @inbounds for k in over_elems_real(grid)
      a_k = q[:a, k, i]
      w_env = q[:w, k, en]
      ρ_k = tmp[:ρ_0, k]
      w_i = q[:w, k, i]
      # dwdz = ∇_z_upwind(q[:w, k-1:k+1, i], q[:w, k-1:k+1, i], grid)
      dwdz = (q[:w, k+1, i]-q[:w, k-1, i])*grid.Δzi*0.5
      B_k = tmp[:buoy, k, i]
      ρa_k = ρ_k * a_k

      nh_press_buoy = - ρa_k * B_k * model.buoy_coef
      nh_pressure_adv = ρa_k * model.nh_avd_coef*w_i*dwdz
      nh_pressure_drag = - ρa_k * model.nh_drag_coef * (w_i - w_env)*abs(w_i - w_env)/max(UpdVar[i].cloud.updraft_top, 500.0)

      tmp[:nh_press, k, i] = nh_press_buoy + nh_pressure_adv + nh_pressure_drag
    end
  end
end

function compute_tke_pressure!(grid::Grid{FT}, UpdVar, q, tmp, tmp_O2, cv, params, model::SCAMPyPressure) where FT
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

function compute_tke_pressure!(grid::Grid{FT}, UpdVar, q, tmp, tmp_O2, cv, params, model::NormalMode) where FT
  gm, en, ud, sd, al = allcombinations(q)
  @inbounds for k in over_elems_real(grid)
    tmp_O2[cv][:press, k] = FT(0)
    @inbounds for i in ud
      a_k = q[:a, k, i]
      w_env = q[:w, k, en]
      ρ_k = tmp[:ρ_0, k]
      w_i = q[:w, k, i]
      dwdz = ∇_z_upwind(q[:w, Cut(k), i], q[:w, Cut(k), i], grid)
      B_k = tmp[:buoy, k, i]
      ρa_k = ρ_k * a_k

      nh_press_buoy = - ρa_k * B_k * model.buoy_coef
      nh_pressure_adv = ρa_k * model.nh_avd_coef*w_i*dwdz
      nh_pressure_drag = - ρa_k * model.nh_drag_coef * (w_i - w_env)*abs(w_i - w_env)/max(UpdVar[i].cloud.updraft_top, 500.0)
      tmp_O2[cv][:press, k] += (w_env-w_i)* nh_pressure_drag
    end
  end
end
