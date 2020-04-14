##### Entrainment-Detrainment models

abstract type EntrDetrModel end
using Printf
struct BOverW2{FT} <: EntrDetrModel
  ε_factor::FT
  δ_factor::FT
end

struct RH_Diff{FT} <: EntrDetrModel
  c_ε::FT
  c_δ::FT
  β_power::FT
  μ_sigmoid::FT
  upd_mixing_frac::FT
  turb_entr_fac::FT
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
      c_ε = model.c_ε
      c_turb = model.turb_entr_fac
      if tmp[:q_liq, k, en]+tmp[:q_liq, k, i]>0.0
        c_δ = model.c_δ
      else
        c_δ = 0.0
      end

      β = model.β_power
      μ_0 = model.μ_sigmoid
      χ = model.upd_mixing_frac

      ts = ActiveThermoState(param_set, q, tmp, k, en)
      RH_en = relative_humidity(ts)
      ts = ActiveThermoState(param_set, q, tmp, k, i)
      RH_up = relative_humidity(ts)
      b_up = tmp[:buoy, k, i]
      b_en = tmp[:buoy, k, en]
      w_up = max(q[:w, k, i],0.1)
      w_en = q[:w, k, en]
      dw = max(w_up - w_en,0.1)
      db = b_up - b_en
      D_ϵ = 1.0/(1.0+exp(-db/dw/μ_0*(χ - q[:a, k, i]/(q[:a, k, i]+q[:a, k, en]))))
      D_δ = 1.0/(1.0+exp( db/dw/μ_0*(χ - q[:a, k, i]/(q[:a, k, i]+q[:a, k, en]))))
      M_δ = ( max((RH_up^β-RH_en^β),0.0) )^(1.0/β)
      M_ϵ = ( max((RH_en^β-RH_up^β),0.0) )^(1.0/β)
      ϵ_dyn = abs(db/dw)/w_up*(c_ε*D_ϵ + c_δ*M_δ)
      δ_dyn = abs(db/dw)/w_up*(c_ε*D_δ + c_δ*M_ϵ)
      ϵ_turb = 2.0#*q[:a, k, i]*c_turb*sqrt(max(q[:tke, k, gm],0.0)) / (q[:w, k, i]*q[:a, k, i]*0.2*UpdVar[i].cloud.updraft_top)
      tmp[:ε_model, k, i] = ϵ_dyn + ϵ_turb
      tmp[:δ_model, k, i] = δ_dyn + ϵ_turb
    end
    tmp[:ε_model, k_1, i] = 2 * Δzi
    tmp[:δ_model, k_1, i] = FT(0)
  end
end

function compute_entrainment_detrainment!(grid::Grid{FT}, UpdVar, tmp, q, params, model::BOverW2) where FT
  gm, en, ud, sd, al = allcombinations(q)
  Δzi = grid.Δzi
  k_1 = first_interior(grid, Zmin())
  @inbounds for i in ud
    zi = UpdVar[i].cloud.base
    @inbounds for k in over_elems_real(grid)
      buoy = tmp[:buoy, k, i]
      w = q[:w, k, i]
      if grid.zc[k] >= zi
        detr_sc = 4.0e-3 + 0.12 *abs(min(buoy,0.0)) / max(w * w, 1e-2)
      else
        detr_sc = FT(0)
      end
      entr_sc = 0.12 * max(buoy, FT(0) ) / max(w * w, 1e-2)
      tmp[:ε_model, k, i] = entr_sc * model.ε_factor
      tmp[:δ_model, k, i] = detr_sc * model.δ_factor
    end
    tmp[:ε_model, k_1, i] = 2 * Δzi
    tmp[:δ_model, k_1, i] = FT(0)
  end
end

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
