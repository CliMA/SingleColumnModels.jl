##### Entrainment-Detrainment models

abstract type EntrDetrModel end
using Printf
# struct BOverW2{FT} <: EntrDetrModel
#   ε_factor::FT
#   δ_factor::FT
# end

struct RH_Diff{FT} <: EntrDetrModel
  ε_factor::FT
  δ_factor::FT
  δ_power::FT
  μ_sigmoid::FT
  upd_mixing_frac::FT
end

"""
    compute_entrainment_detrainment!

Define entrainment and detrainment fields

 - `tmp[:ε_model, k, i]`
 - `tmp[:δ_model, k, i]`

 for all `k` and all `i`
"""
function compute_entrainment_detrainment! end

function compute_entrainment_detrainment!(grid::Grid{FT}, UpdVar, tmp, q, params, model::RH_Diff) where FT
  gm, en, ud, sd, al = allcombinations(q)
  Δzi = grid.Δzi
  k_1 = first_interior(grid, Zmin())
  @inbounds for i in ud
    @inbounds for k in over_elems_real(grid)
      @unpack params param_set
      c_ent = model.ε_factor
      if tmp[:q_liq, k, en]+tmp[:q_liq, k, i]>0.0
        c_det = model.δ_factor
      else
        c_det = 0.0
      end
      if tmp[:q_liq, k, en]+tmp[:q_liq, k, i]>0.0
        c_det = model.δ_factor
      else
        c_det = 0.0
      end

      β = model.δ_power
      μ = model.μ_sigmoid
      χ = model.upd_mixing_frac

      ts = ActiveThermoState(param_set, q, tmp, k, en)
      RH_en = relative_humidity(ts)
      ts = ActiveThermoState(param_set, q, tmp, k, i)
      RH_up = relative_humidity(ts)
      b_up = tmp[:buoy, k, i]
      b_en = tmp[:buoy, k, en]
      w_up = max(q[:w, k, i],0.001)
      w_en = q[:w, k, en]
      dw = max(w_up - w_en,0.001)
      db = b_up - b_en
      logistic_e = 1.0/(1.0+exp(-μ*db/dw*(χ - q[:a, k, i]/(q[:a, k, i]+q[:a, k, en]))))
      logistic_d = 1.0/(1.0+exp( μ*db/dw*(χ - q[:a, k, i]/(q[:a, k, i]+q[:a, k, en]))))
      moisture_deficit_d = (max(RH_up^β-RH_en^β,0.0))^(1.0/β)
      moisture_deficit_e = (max(RH_en^β-RH_up^β,0.0))^(1.0/β)
      tmp[:ε_model, k, i] = abs(db/dw)/w_up*(c_ent*logistic_e + c_det*moisture_deficit_e)
      tmp[:δ_model, k, i] = abs(db/dw)/w_up*(c_ent*logistic_d + c_det*moisture_deficit_e)
      @printf("%F %f %f %f %f\n", tmp[:ε_model, k, i], tmp[:δ_model, k, i], q[:a, k, i], w_up, b_up)
    end
    tmp[:ε_model, k_1, i] = 2 * Δzi
    tmp[:δ_model, k_1, i] = FT(0)
  end
end

# function compute_entrainment_detrainment!(grid::Grid{FT}, UpdVar, tmp, q, params, model::BOverW2) where FT
#   gm, en, ud, sd, al = allcombinations(q)
#   Δzi = grid.Δzi
#   k_1 = first_interior(grid, Zmin())
#   @inbounds for i in ud
#     zi = UpdVar[i].cloud.base
#     @inbounds for k in over_elems_real(grid)
#       buoy = tmp[:buoy, k, i]
#       w = q[:w, k, i]
#       if grid.zc[k] >= zi
#         detr_sc = 4.0e-3 + 0.12 *abs(min(buoy,0.0)) / max(w * w, 1e-2)
#       else
#         detr_sc = FT(0)
#       end
#       entr_sc = 0.12 * max(buoy, FT(0) ) / max(w * w, 1e-2)
#       tmp[:ε_model, k, i] = entr_sc * model.ε_factor
#       tmp[:δ_model, k, i] = detr_sc * model.δ_factor
#     end
#     tmp[:ε_model, k_1, i] = 2 * Δzi
#     tmp[:δ_model, k_1, i] = FT(0)
#   end
# end

function compute_cv_entr!(grid::Grid{FT}, q, tmp, tmp_O2, ϕ, ψ, cv, tke_factor) where FT
  gm, en, ud, sd, al = allcombinations(q)
  @inbounds for k in over_elems_real(grid)
    tmp_O2[cv][:entr_gain, k] = FT(0)
    @inbounds for i in ud
      Δϕ = q[ϕ, k, i] - q[ϕ, k, en]
      Δψ = q[ψ, k, i] - q[ψ, k, en]
      tmp_O2[cv][:entr_gain, k] += tke_factor*q[:a, k, i] * abs(q[:w, k, i]) * tmp[:δ_model, k, i] * Δϕ * Δψ
    end
    tmp_O2[cv][:entr_gain, k] *= tmp[:ρ_0, k]
  end
end
