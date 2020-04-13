##### Entrainment-Detrainment models

abstract type EntrDetrModel end

struct BOverW2{FT} <: EntrDetrModel
  ε_factor::FT
  δ_factor::FT
  δ_power::FT
end

"""
    compute_entrainment_detrainment!

Define entrainment and detrainment fields

 - `tmp[:ε_model, k, i]`
 - `tmp[:δ_model, k, i]`

 for all `k` and all `i`
"""
function compute_entrainment_detrainment! end

function compute_entrainment_detrainment!(grid::Grid{FT}, UpdVar, tmp, q, params, model::BOverW2) where FT
  gm, en, ud, sd, al = allcombinations(q)
  Δzi = grid.Δzi
  k_1 = first_interior(grid, Zmin())
  beta = model.δ_power
  c_ent = model.ε_factor
  c_det = model.δ_factor
  @inbounds for k in over_elems_real(grid)
    RH_en = relative_humidity(q, tmp, k, en)
    w_env = q[:w, k, en]
    b_env = tmp[:buoy, k, env]
    @inbounds for i in ud
      db = tmp[:buoy, k, i] - b_env
      w_up = q[:w, k, i]
      dw = w_up - w_env
      RH_up = relative_humidity(q, tmp, k, i)
      if dw < 0.0
        dw = dw -  0.001
      else
        dw = dw + 0.001
      end
      if ([:ql, k, i]+[:ql, k, en])==0.0
        c_det = 0.0
      else
        c_det = 1.0
      end
      moisture_deficit_d = (max(RH_upd^2.0-RH_env^2.0,0.0))^(1.0/2.0)
      moisture_deficit_e = (max(RH_env^2.0-RH_upd^2.0,0.0))^(1.0/2.0)
      tmp[:ε_model, k, i] = abs(db/dw)/max(w_up*w_up,1e-4)*(c_ent*logistic_e+c_det*moisture_deficit_e)
      tmp[:δ_model, k, i] = abs(db/dw)/max(w_up*w_up,1e-4)*(c_ent*logistic_d+c_det*moisture_deficit_e)
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
