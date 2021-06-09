##### Buoyancy models

@inline buoyancy(param_set, α_0::FT, α::FT) where {FT} =
    FT(grav(param_set)) * (α - α_0) / α_0

"""
    compute_buoyancy!

Define buoyancy field

 - `aux[:buoy, k, i]`

 for all `k` and all `i`
"""
function compute_buoyancy! end

function compute_buoyancy!(grid, q, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set = params
    @inbounds for i in (ud..., en)
        @inbounds for k in over_elems_real(grid)
            q_tot = q[:q_tot, k, i]
            q_liq = aux[:q_liq, k, i]
            T = aux[:T, k, i]
            α_i = specific_volume(
                param_set,
                T,
                aux[:p_0, k],
                PhasePartition(q_tot, q_liq),
            )
            aux[:buoy, k, i] = buoyancy(param_set, aux[:α_0, k], α_i)
        end
    end

    # Filter buoyancy
    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            weight = aux[:HVSD_a, k, i]
            aux[:buoy, k, i] =
                weight * aux[:buoy, k, i] + (1 - weight) * aux[:buoy, k, en]
        end
    end

    # Subtract grid mean buoyancy
    @inbounds for k in over_elems_real(grid)
        aux[:buoy, k, gm] = sum([q[:a, k, i] * aux[:buoy, k, i] for i in sd])
        @inbounds for i in sd
            aux[:buoy, k, i] -= aux[:buoy, k, gm]
        end
    end
end

function compute_tke_buoy!(
    grid::Grid{FT},
    q,
    aux,
    aux_O2,
    cv,
    params,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set = params

    _R_v::FT = R_v(param_set)
    _R_d::FT = R_d(param_set)
    _grav::FT = grav(param_set)

    # Note that source terms at the first interior point are not really used because that is where tke boundary condition is
    # enforced (according to MO similarity). Thus here I am being sloppy about lowest grid point
    @inbounds for k in over_elems_real(grid)
        q_tot_dry = aux[:q_tot_dry, k]
        θ_dry = aux[:θ_dry, k]
        t_cloudy = aux[:t_cloudy, k]
        q_vap_cloudy = aux[:q_vap_cloudy, k]
        q_tot_cloudy = aux[:q_tot_cloudy, k]
        θ_cloudy = aux[:θ_cloudy, k]
        p_0 = aux[:p_0, k]

        lh = latent_heat_vapor(param_set, t_cloudy)
        cpm = cp_m(param_set, PhasePartition(q_tot_cloudy))
        grad_θ_liq = ∇_dn(q[:θ_liq, Cut(k), en], grid)
        grad_q_tot = ∇_dn(q[:q_tot, Cut(k), en], grid)

        prefactor = _R_d * TD.exner_given_pressure(param_set, p_0) / p_0
        ε_vi = _R_v / _R_d

        d_alpha_θ_liq_dry = prefactor * (1 + (ε_vi - 1) * q_tot_dry)
        d_alpha_q_tot_dry = prefactor * θ_dry * (ε_vi - 1)
        CF_env = aux[:CF, k]

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

        K_h_k = aux[:K_h, k, gm]
        term_1 = -K_h_k * grad_θ_liq * d_alpha_θ_liq_total
        term_2 = -K_h_k * grad_q_tot * d_alpha_q_tot_total

        # TODO - check
        aux_O2[cv][:buoy, k] =
            _grav / aux[:α_0, k] *
            q[:a, k, en] *
            aux[:ρ_0, k] *
            (term_1 + term_2)
    end
end
