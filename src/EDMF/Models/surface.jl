##### Surface models

abstract type SurfaceModel end

struct SurfaceFixedFlux{FT} <: SurfaceModel
    T::FT
    P::FT
    q_tot::FT
    shf::FT
    lhf::FT
    Tsurface::FT
    ρ_0_surf::FT
    α_0_surf::FT
    ρq_tot_flux::FT
    ρθ_liq_flux::FT
    ustar::FT
    windspeed_min::FT
    tke_tol::FT
    area::FT
end

function SurfaceFixedFlux(
    param_set;
    T::FT,
    P,
    q_tot,
    ustar,
    windspeed_min,
    tke_tol,
    area,
) where {FT}
    q_pt = PhasePartition(q_tot)
    ρ_0_surf = air_density(param_set, T, P, q_pt)
    α_0_surf = 1 / ρ_0_surf
    Tsurface = 299.1 * TD.exner_given_pressure(param_set, P)
    lhf = 5.2e-5 * ρ_0_surf * latent_heat_vapor(param_set, Tsurface)
    shf = 8.0e-3 * cp_m(param_set, q_pt) * ρ_0_surf
    ρ_tflux = shf / cp_m(param_set, q_pt)
    ρq_tot_flux = lhf / latent_heat_vapor(param_set, Tsurface)
    ρθ_liq_flux = ρ_tflux / TD.exner_given_pressure(param_set, P)
    return SurfaceFixedFlux{FT}(
        T,
        P,
        q_tot,
        shf,
        lhf,
        Tsurface,
        ρ_0_surf,
        α_0_surf,
        ρq_tot_flux,
        ρθ_liq_flux,
        ustar,
        windspeed_min,
        tke_tol,
        area,
    )
end

function SurfaceFixedFlux_dycoms(
    param_set;
    T::FT,
    shf,
    lhf,
    Tsurface,
    qsurface,
    P,
    q_tot,
    ustar,
    windspeed_min,
    tke_tol,
    area,
) where {FT}
    q_pt = PhasePartition(q_tot)
    ρ_0_surf = air_density(param_set, T, P, q_pt)
    α_0_surf = 1 / ρ_0_surf
    ρ_tflux = shf / cp_m(param_set, q_pt)
    ρq_tot_flux = lhf / latent_heat_vapor(param_set, Tsurface)
    ρθ_liq_flux = shf / cp_m(param_set, q_pt) / ρ_0_surf
    θ_surface = Tsurface / TD.exner_given_pressure(param_set, P)

    num =
        ρθ_liq_flux +
        (eps_vi - 1.0) *
        (θ_surface * ρq_tot_flux / ρ_0_surf + qsurface * ρθ_liq_flux)
    denom = (θ_surface * (1.0 + (eps_vi - 1) * qsurface))
    bflux = FT(grav(param_set)) * num / denom

    return SurfaceFixedFlux{FT}(
        T,
        P,
        q_tot,
        shf,
        lhf,
        Tsurface,
        ρ_0_surf,
        α_0_surf,
        ρq_tot_flux,
        ρθ_liq_flux,
        ustar,
        windspeed_min,
        tke_tol,
        area,
    )
end
