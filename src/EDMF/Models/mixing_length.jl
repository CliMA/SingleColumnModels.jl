##### Mixing length models

abstract type MixingLengthModel end

#####
##### Types
#####

"""
    ConstantMixingLength

TODO: add docs
"""
struct ConstantMixingLength{FT} <: MixingLengthModel
  value::FT
end

"""
    SCAMPyMixingLength

TODO: add docs
"""
struct SCAMPyMixingLength{FT} <: MixingLengthModel
  a_L::StabilityDependentParam{FT}
  b_L::StabilityDependentParam{FT}
end

"""
    IgnaciosMixingLength

TODO: add docs
"""
struct IgnaciosMixingLength{FT} <: MixingLengthModel
  a_L::StabilityDependentParam{FT}
  b_L::StabilityDependentParam{FT}
  c_K::FT
  c_ε::FT
  c_w::FT
  ω_1::FT
  ω_2::FT
  Prandtl_neutral::FT
  function IgnaciosMixingLength(a_L::StabilityDependentParam{FT},
                                b_L::StabilityDependentParam{FT},
                                c_K::FT,
                                c_ε::FT,
                                c_w::FT,
                                ω_1::FT,
                                Prandtl_neutral::FT) where FT
    return new{FT}(a_L, b_L, c_K, c_ε, c_w, ω_1, ω_1 + 1, Prandtl_neutral)
  end
end

######
###### Helper funcs
######

"""
    compute_mixing_τ

TODO: add docs
"""
compute_mixing_τ(zi::FT, wstar::FT) where FT = zi / (wstar + FT(0.001))

"""
    ϕ_m

TODO: add docs
"""
ϕ_m(ξ, a_L, b_L) = (1+a_L*ξ)^(-b_L)

######
###### Compute mixing length methods
######

"""
    compute_mixing_length!

Define mixing length

 - `tmp[:l_mix, k, i]`

for all `k` and all `i`
"""
function compute_mixing_length! end

function compute_mixing_length!(grid::Grid{FT}, q, tmp, params, model::ConstantMixingLength) where FT
  gm, en, ud, sd, al = allcombinations(q)
  @inbounds for k in over_elems(grid)
    tmp[:l_mix, k, gm] = model.value
  end
end

function compute_mixing_length!(grid::Grid{FT}, q, tmp, params, model::SCAMPyMixingLength{FT}) where FT
  @unpack params obukhov_length zi wstar k_Karman
  gm, en, ud, sd, al = allcombinations(q)
  τ = compute_mixing_τ(zi, wstar)
  a_L = model.a_L(obukhov_length)
  b_L = model.b_L(obukhov_length)
  @inbounds for k in over_elems_real(grid)
    l1 = τ * sqrt(max(q[:tke, k, en], FT(0)))
    z = grid.zc[k]
    ξ = z/obukhov_length
    l2 = k_Karman * z * ϕ_m(ξ, a_L, b_L)
    tmp[:l_mix, k, gm] = max( 1/(1/max(l1,1e-10) + 1/l2), FT(1e-3))
  end
end

function compute_mixing_length!(grid::Grid{FT}, q, tmp, params, model::IgnaciosMixingLength{FT}) where FT
  @unpack params obukhov_length SurfaceModel k_Karman param_set
  gm, en, ud, sd, al = allcombinations(q)
  ustar = SurfaceModel.ustar
  L = Vector(undef, 3)
  a_L = model.a_L(obukhov_length)
  b_L = model.b_L(obukhov_length)
  @inbounds for k in over_elems_real(grid)
    TKE_k = max(q[:tke, k, en], FT(0))
    θ_ρ = tmp[:θ_ρ, k, gm]
    z = grid.zc[k]
    ts_dual = ActiveThermoState(param_set, q, tmp, Dual(k), gm)
    θ_ρ_dual = virtual_pottemp.(ts_dual)
    ∇θ_ρ = ∇_z_flux(θ_ρ_dual, grid)
    buoyancy_freq = FT(grav(param_set))*∇θ_ρ/θ_ρ
    L[1] = sqrt(model.c_w*TKE_k)/buoyancy_freq
    ξ = z/obukhov_length
    κ_star = ustar/sqrt(TKE_k)
    L[2] = k_Karman*z/(model.c_K*κ_star*ϕ_m(ξ, a_L, b_L))
    S_squared = ∇_z_flux(q[:u, Dual(k), gm], grid)^2 +
                ∇_z_flux(q[:v, Dual(k), gm], grid)^2 +
                ∇_z_flux(q[:w, Dual(k), en], grid)^2
    R_g = tmp[:∇buoyancy, k, gm]/S_squared
    if unstable(obukhov_length)
      Pr_z = model.Prandtl_neutral
    else
      discriminant1 = -4*R_g+(1+model.ω_2*R_g)^2
      Pr_z = model.Prandtl_neutral*(1+model.ω_2*R_g - sqrt(discriminant1))/(2*R_g)
    end
    discriminant2 = S_squared - tmp[:∇buoyancy, k, gm]/Pr_z
    L[3] = sqrt(model.c_ε/model.c_K)*sqrt(TKE_k)*1/sqrt(max(discriminant2, 1e-2))
    tmp[:l_mix, k, gm] = sum([L[j]*exp(-L[j]) for j in 1:3])/sum([exp(-L[j]) for j in 1:3])
  end
end
