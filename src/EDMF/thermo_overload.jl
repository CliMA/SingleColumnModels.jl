using RootSolvers

"""
    specific_volume(T, p[, q::PhasePartition])
The (moist-)air specific volume from the equation of
state (ideal gas law) where
 - `T` air temperature
 - `p` pressure
and, optionally,
 - `q` [`PhasePartition`](@ref). Without this argument the results are that of dry air.
Old method (used in SCAMPY) for specific volume
TODO: Remove this method / synchronize with `specific_volume(ts::ThermodynamicState)`
"""
specific_volume(param_set::PS, T::FT, p::FT, q::PhasePartition{FT}=PhasePartition{FT}(FT(0), FT(0), FT(0))) where {FT<:Real, PS} =
  (gas_constant_air(param_set, q) * T) / p

"""
    buoyancy_flux(shf::FT, lhf::FT, T_b::FT, q_tot::FT, α_0::FT) where {FT<:Real}
Old method (used in SCAMPY) for buoyancy flux
TODO: Remove this method.
"""
function buoyancy_flux(param_set::PS, shf::FT, lhf::FT, T_b::FT, q_tot::FT, α_0::FT) where {FT<:Real, PS}
  cp_ = cp_m(param_set, PhasePartition(q_tot))
  lv = latent_heat_vapor(param_set, T_b)
  return (FT(grav(param_set)) * α_0 / cp_ / T_b * (shf + ((FT(R_v(param_set)) / FT(R_d(param_set)))-1) * cp_ * T_b * lhf /lv))
end

"""
    air_temperature_from_liquid_ice_pottemp_old(θ_liq_ice::FT, p::FT, q::PhasePartition{FT}=q_pt_0(FT))
Old method (used in SCAMPY) for air temperature from liquid ice potential temperature
TODO: remove this method / synchronize with `air_temperature_from_liquid_ice_pottemp`
"""
air_temperature_from_liquid_ice_pottemp_old(param_set::PS, θ_liq_ice::FT, p::FT, q::PhasePartition{FT}=Thermodynamics.q_pt_0(FT)) where {FT<:Real, PS} =
  θ_liq_ice*TD.exner_given_pressure(param_set, p, q) + (FT(LH_v0(param_set))*q.liq + FT(LH_s0(param_set))*q.ice) / cp_m(param_set, q)

air_temperature_from_liquid_ice_pottemp_old(param_set::PS, θ_liq_ice::FT, p::FT) where {FT<:Real, PS} =
  air_temperature_from_liquid_ice_pottemp_old(param_set, θ_liq_ice, p, Thermodynamics.q_pt_0(FT))

"""
    saturation_adjustment_q_tot_θ_liq_ice(θ_liq_ice, q_tot, ρ, p)
Compute the temperature that is consistent with
 - `θ_liq_ice` liquid-ice potential temperature
 - `q_tot` total specific humidity
 - `ρ` density
 - `p` pressure
See also [`saturation_adjustment`](@ref).
Old method (used in SCAMPY) for `saturation_adjustment_q_tot_θ_liq_ice`
TODO: remove this method / synchronize with `saturation_adjustment_q_tot_θ_liq_ice`
"""
function saturation_adjustment_q_tot_θ_liq_ice_old(param_set::PS,
                                                   θ_liq_ice::FT,
                                                   q_tot::FT,
                                                   ρ::FT,
                                                   p::FT) where {FT<:Real, PS}
  T_1 = air_temperature_from_liquid_ice_pottemp_old(param_set, θ_liq_ice, p) # Assume all vapor
  q_v_sat = TD.q_vap_saturation(param_set, T_1, ρ, PhaseEquil)
  T_sol = T_1
  if q_tot <= q_v_sat # If not saturated
  else  # If saturated, iterate
    q_pt_ice = PhasePartition(q_tot, FT(0), q_tot) # Assume all ice
    T_2 = air_temperature_from_liquid_ice_pottemp_old(param_set, θ_liq_ice, p, q_pt_ice)
    function eos(T)
      q_vap_sat = TD.q_vap_saturation(param_set, T, ρ, PhaseEquil)
      q_pt = PhasePartition(q_vap_sat)
      θ_dry = TD.dry_pottemp_given_pressure(param_set, T, p, q_pt)
      exp_fact = exp(-latent_heat_vapor(param_set, T)/(T*cp_d(param_set))*(q_tot - q_vap_sat)/(1.0-q_tot))
      return θ_liq_ice - θ_dry * exp_fact
    end

    # TODO/TOFIX: SCAMPY implementation
    tol, maxiter = ResidualTolerance(FT(1e-3)), 10
    sol = find_zero(T -> eos(T), SecantMethod(T_1, T_2), CompactSolution(), tol, maxiter)
    T_sol = sol.root
  end
  return T_sol
end

"""
    LiquidIcePotTempSHumEquil(θ_liq_ice, q_tot, ρ, p)
Constructs a [`PhaseEquil`](@ref) thermodynamic state from:
 - `θ_liq_ice` - liquid-ice potential temperature
 - `q_tot` - total specific humidity
 - `ρ` - density
 - `p` - pressure
"""
function LiquidIcePotTempSHumEquil_old(param_set::PS, θ_liq_ice::FT, q_tot::FT, ρ::FT, p::FT) where {FT<:Real,PS}
    T = saturation_adjustment_q_tot_θ_liq_ice_old(param_set, θ_liq_ice, q_tot, ρ, p)
    q = PhasePartition_equil(param_set, T, ρ, q_tot, PhaseEquil)
    e_int = internal_energy(param_set, T, q)
    return PhaseEquil{FT,PS}(param_set, e_int, ρ, q_tot, T)
end
