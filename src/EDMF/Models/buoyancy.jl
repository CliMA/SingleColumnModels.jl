##### Buoyancy models

@inline buoyancy(param_set, α_0::FT, α::FT) where {FT} =
    FT(grav(param_set)) * (α - α_0) / α_0

"""
    compute_buoyancy!

Define buoyancy field

 - `tmp[:buoy, k, i]`

 for all `k` and all `i`
"""
function compute_buoyancy! end

function compute_buoyancy!(grid, q, tmp, params)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack params param_set
    @inbounds for i in (ud..., en)
        @inbounds for k in over_elems_real(grid)
            q_tot = q[:q_tot, k, i]
            q_liq = tmp[:q_liq, k, i]
            T = tmp[:T, k, i]
            α_i = specific_volume(
                param_set,
                T,
                tmp[:p_0, k],
                PhasePartition(q_tot, q_liq),
            )
            tmp[:buoy, k, i] = buoyancy(param_set, tmp[:α_0, k], α_i)
        end
    end

    # Filter buoyancy
    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            weight = tmp[:HVSD_a, k, i]
            tmp[:buoy, k, i] =
                weight * tmp[:buoy, k, i] + (1 - weight) * tmp[:buoy, k, en]
        end
    end

    # Subtract grid mean buoyancy
    @inbounds for k in over_elems_real(grid)
        tmp[:buoy, k, gm] = sum([q[:a, k, i] * tmp[:buoy, k, i] for i in sd])
        @inbounds for i in sd
            tmp[:buoy, k, i] -= tmp[:buoy, k, gm]
        end
    end
end

function compute_tke_buoy!(
    grid::Grid{FT},
    q,
    tmp,
    tmp_O2,
    cv,
    params,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @unpack params param_set

    _R_v::FT = R_v(param_set)
    _R_d::FT = R_d(param_set)
    _grav::FT = grav(param_set)

    # Note that source terms at the first interior point are not really used because that is where tke boundary condition is
    # enforced (according to MO similarity). Thus here I am being sloppy about lowest grid point
    @inbounds for k in over_elems_real(grid)
        q_tot_dry = tmp[:q_tot_dry, k]
        θ_dry = tmp[:θ_dry, k]
        t_cloudy = tmp[:t_cloudy, k]
        q_vap_cloudy = tmp[:q_vap_cloudy, k]
        q_tot_cloudy = tmp[:q_tot_cloudy, k]
        θ_cloudy = tmp[:θ_cloudy, k]
        p_0 = tmp[:p_0, k]

        lh = latent_heat_vapor(param_set, t_cloudy)
        cpm = cp_m(param_set, PhasePartition(q_tot_cloudy))
        grad_θ_liq = ∇_dn(q[:θ_liq, Cut(k), en], grid)
        grad_q_tot = ∇_dn(q[:q_tot, Cut(k), en], grid)

        prefactor = _R_d * TD.exner_given_pressure(param_set, p_0) / p_0
        ε_vi = _R_v / _R_d

        d_alpha_θ_liq_dry = prefactor * (1 + (ε_vi - 1) * q_tot_dry)
        d_alpha_q_tot_dry = prefactor * θ_dry * (ε_vi - 1)
        CF_env = tmp[:CF, k]

        if CF_env > 0
            d_alpha_θ_liq_cloudy = (
                prefactor * (
                    1 + ε_vi * (1 + lh / _R_v / t_cloudy) * q_vap_cloudy -
                    q_tot_cloudy
                ) / (
                    1 +
                    lh * lh / cpm / _R_v / t_cloudy / t_cloudy * q_vap_cloudy
                )
            )
            d_alpha_q_tot_cloudy =
                (lh / cpm / t_cloudy * d_alpha_θ_liq_cloudy - prefactor) *
                θ_cloudy
        else
            d_alpha_θ_liq_cloudy = 0
            d_alpha_q_tot_cloudy = 0
        end

        d_alpha_θ_liq_total =
            (CF_env * d_alpha_θ_liq_cloudy + (1 - CF_env) * d_alpha_θ_liq_dry)
        d_alpha_q_tot_total =
            (CF_env * d_alpha_q_tot_cloudy + (1 - CF_env) * d_alpha_q_tot_dry)

        K_h_k = tmp[:K_h, k, gm]
        term_1 = -K_h_k * grad_θ_liq * d_alpha_θ_liq_total
        term_2 = -K_h_k * grad_q_tot * d_alpha_q_tot_total

        # TODO - check
        tmp_O2[cv][:buoy, k] =
            _grav / tmp[:α_0, k] *
            q[:a, k, en] *
            tmp[:ρ_0, k] *
            (term_1 + term_2)
    end
end
