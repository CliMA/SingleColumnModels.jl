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

struct MinimumDissipation{FT} <: MixingLengthModel
  a_L::StabilityDependentParam{FT}
  b_L::StabilityDependentParam{FT}
  c_K::FT
  c_ε::FT
  c_w::FT
  ω_1::FT
  ω_2::FT
  Prandtl_neutral::FT
  function MinimumDissipation(a_L::StabilityDependentParam{FT},
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
    # ϕ_m(ξ, a_L, b_L)) = (1+a_L*ξ)^(-b_L)
    # L[2] = k_Karman*z/(model.c_K*κ_star*(1+a_L*ξ)^(-b_L)
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

function compute_mixing_length!(grid::Grid{FT}, q, tmp, params, model::MinimumDissipation{FT}) where FT
  @unpack params obukhov_length SurfaceModel k_Karman param_set
  gm, en, ud, sd, al = allcombinations(q)
  ustar = SurfaceModel.ustar
  L = Vector(undef, 3)
  a_L = model.a_L(obukhov_length)
  b_L = model.b_L(obukhov_length)
  @inbounds for k in over_elems_real(grid)
    z = grid.zc[k]
    if obukhov_length < 0.0 #unstable
      L[2] = k_Karman * z/(sqrt(max(q[:tke, gw, en], FT(0))/ustar/ustar)* model.c_K) * fmin(
         (1 - 100 * z/obukhov_length)^0.2, 1/k_Karman ) # m ake sure it TKE in first interior
    else # neutral or stable
      L[2] = k_Karman * z/(sqrt(max(q[:tke, gw, en], FT(0))/ustar/ustar)*model.c_K)
    end
    TKE_k = max(q[:tke, k, en], FT(0))

    θ_ρ = tmp[:θ_ρ, k, gm]
    z = grid.zc[k]
    ts_dual = ActiveThermoState(param_set, q, tmp, Dual(k), gm)
    θ_ρ_dual = virtual_pottemp.(ts_dual)
    ∇θ_ρ = ∇_z_flux(θ_ρ_dual, grid)
    buoyancy_freq = FT(grav(param_set))*∇θ_ρ/θ_ρ
    if buoyancy_freq>0.0
      L[1] = sqrt(model.c_w*TKE_k)/buoyancy_freq # check if this can be singular
    else
      L[1] = 1e-6
    end
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
      Pr_z[k] = model.Prandtl_neutral*(2*R_g/
                        (1+(53/13)*R_g -sqrt( (1+(53/130)*R_g)^2 - 4*R_g ) ) )
    end

    # missing grad_b_thl
    grad_b_thl = 1

    # Production/destruction terms
    a = model.c_ε*(S_squared - grad_b_thl/Pr_z[k] - grad_b_qt/Pr_z[k])* sqrt(TKE_k)
    # Dissipation term
    c_neg = model.c_ε*TKE_k*sqrt(TKE_k)
    b[k] = 0.0
    @inbounds for i in ud
        b[k] += q[:a, k, i]*q[:w, k, i]*tmp[:δ_model, k, i]/q[:a, k, en]*((q[:w, k, i]-q[:w, k, en])^2/2.0-TKE_k) - q[:a, k, i]*q[:w, k, i]*(
            q[:w, k, i]-q[:w, k, en])*tmp[:εt_model, k, i]*q[:w, k, en]/q[:a, k, en]
    end

    if abs(a) > m_eps && 4*a*c_neg > - b[k]^2
              l_entdet[k] = max( -b[k]/2.0/a + sqrt(b[k]^2 + 4*a*c_neg)/2/a, 0)
    elseif abs(a) < m_eps && abs(b[k]) > m_eps
              l_entdet[k] = c_neg/b[k]
    end
    L[3] = l_entdet[k]
    tmp[:l_mix, k, gm] = sum([L[j]*exp(-L[j]) for j in 1:3])/sum([exp(-L[j]) for j in 1:3])
  end
end