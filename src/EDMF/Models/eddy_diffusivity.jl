##### Eddy-diffusivity models

abstract type EddyDiffusivityModel end

Base.@kwdef struct SCAMPyEddyDiffusivity{FT} <: EddyDiffusivityModel
  tke_ed_coeff::FT=FT(0)
  prandtl_number::FT=FT(0)
end

"""
    compute_eddy_diffusivities_tke!

Define eddy-diffusivity fields

 - `tmp[:K_m, k, i]`
 - `tmp[:K_h, k, i]`

 for all `k` and all `i`
"""
function compute_eddy_diffusivities_tke! end

function compute_eddy_diffusivities_tke!(grid::Grid{FT}, q, tmp, params, model::SCAMPyEddyDiffusivity) where FT
  gm, en, ud, sd, al = allcombinations(q)
  @inbounds for k in over_elems_real(grid)
    l_mix = tmp[:l_mix, k, gm]
    K_m_k = model.tke_ed_coeff * l_mix * sqrt(max(q[:tke, k, en], FT(0)))
    tmp[:K_m, k, gm] = K_m_k
    tmp[:K_h, k, gm] = K_m_k / model.prandtl_number
  end
end
