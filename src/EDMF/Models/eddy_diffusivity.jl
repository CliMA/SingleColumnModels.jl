##### Eddy-diffusivity models

abstract type EddyDiffusivityModel end
struct SCAMPyEddyDiffusivity{FT} <: EddyDiffusivityModel
  tke_ed_coeff::FT
end

function compute_eddy_diffusivities_tke!(grid::Grid{FT}, q, tmp, params, model::SCAMPyEddyDiffusivity) where FT
  gm, en, ud, sd, al = allcombinations(q)
  @inbounds for k in over_elems_real(grid)
    l_mix = tmp[:l_mix, k, gm]
    K_m_k = model.tke_ed_coeff * l_mix * sqrt(max(q[:tke, k, en], FT(0)))
    tmp[:K_m, k, gm] = K_m_k
    tmp[:K_h, k, gm] = K_m_k / params[:prandtl_number]
  end
end
