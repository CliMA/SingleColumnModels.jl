##### SGS quadrature
using FastGaussQuadrature
abstract type SubdomainStatistics end

struct SubdomainMean{IT <: Int} <: SubdomainStatistics
    order::IT
end


struct GaussianQuadrature{IT <: Int} <: SubdomainStatistics
    order::IT
end

struct LogNormalQuadrature{IT <: Int} <: SubdomainStatistics
    order::IT
end

"""
    compute_subdomain_statistics!

at the moment only in the en subdomain
"""

function compute_subdomain_statistics! end

function compute_subdomain_statistics!(
    grid::Grid{FT},
    i,
    q,
    aux,
    params,
    model::SubdomainMean,
) where {FT, IT}
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set = params
    @inbounds for k in over_elems_real(grid)
        q_tot = q[:q_tot, k, en]
        ts = ActiveThermoState(param_set, q, aux, k, en)
        T = air_temperature(ts)
        q_liq = PhasePartition(ts).liq
        q_vap = q_tot - q_liq
        θ = dry_pottemp(ts)
        if q_liq > 0
            aux[:CF, k] = FT(1)
            aux[:θ_cloudy, k] = θ
            aux[:t_cloudy, k] = T
            aux[:q_tot_cloudy, k] = q_tot
            aux[:q_vap_cloudy, k] = q_vap
        else
            aux[:CF, k] = FT(0)
            aux[:θ_dry, k] = θ
            aux[:q_tot_dry, k] = q_tot
        end
    end
end

function compute_subdomain_statistics!(
    grid::Grid{FT},
    i,
    q,
    aux,
    params,
    model::GaussianQuadrature,
) where {FT, IT}
    gm, en, ud, sd, al = allcombinations(q)
    abscissas, weights = gausshermite(model.order)
    sqrt2 = sqrt(2.0)
    sqpi_inv = 1 / sqrt(pi)
    env_len = 10
    src_len = 6
    @unpack param_set = params

    inner_env = zeros(env_len)
    outer_env = zeros(env_len)
    inner_src = zeros(src_len)
    outer_src = zeros(src_len)
    i_Sqt = collect(1:src_len)
    i_SH = collect(1:src_len)
    i_Sqt_H = collect(1:src_len)
    i_Sqt_qt = collect(1:src_len)
    i_SH_H = collect(1:src_len)
    i_SH_qt = collect(1:src_len)
    i_ql = collect(1:env_len)
    i_T = collect(1:env_len)
    i_θl = collect(1:env_len)
    i_ρ = collect(1:env_len)
    i_cf = collect(1:env_len)
    i_qt_cld = collect(1:env_len)
    i_qt_dry = collect(1:env_len)
    i_T_cld = collect(1:env_len)
    i_T_dry = collect(1:env_len)
    i_rf = collect(1:env_len)

    @inbounds for k in over_elems_real(grid)
        if (
            q[:cv_q_tot, k, i] > eps(FT) &&
            q[:cv_θ_liq, k, i] > eps(FT) &&
            fabs(q[:cv_θ_liq_q_tot, k, i]) > eps(FT) &&
            q[:q_tot, k, i] > eps(FT) &&
            sqrt(q[:cv_q_tot, k, i]) < q[:q_tot, k, i]
        )
            σ_q = sqrt(q[:cv_q_tot, k, i])
            σ_h = sqrt(q[:cv_θ_liq, k, i])
            corr = max(
                min(q[:cv_θ_liq_q_tot, k, i] / max(σ_h * σ_q, 1e-13), 1),
                -1,
            )
            # limit σ_q to prevent negative qt_hat
            σ_qt_lim = (1e-10 - q[:q_tot, k, i]) / (sqrt2 * abscissas[0])
            # walking backwards to assure your q_t will not be smaller than 1e-10
            # TODO - check
            # TODO - change 1e-13 and 1e-10 to some epislon
            σ_q = min(σ_q, σ_qt_lim)
            qt_covar = σ_q * σ_q
            σ_θ_star = sqrt(max(1 - corr * corr, 0.0)) * σ_h

            # clean outer vectors
            outer_src = 0.0 * outer_src
            outer_env = 0.0 * outer_env

            for j_qt in 1:(model.order)
                qt_hat = q[:q_tot, k, i] + sqrt2 * σ_q * abscissas[j_qt]
                μ_θ_star =
                    q[:θ_lid, k, i] + sqrt2 * corr * σ_h * abscissas[j_qt]
                # clean innner vectors
                inner_src = 0.0 * inner_src
                inner_env = 0.0 * inner_env
                for j_θ in 1:(model.order)
                    θ_hat = sqrt2 * σ_θ_star * abscissas[j_θ] + μ_θ_star
                    ts = DummyThermoState(
                        param_set;
                        θ = θ_hat,
                        q_tot = qt_hat,
                        ρ_0 = aux[:ρ_0, k, i],
                        p_0 = aux[:p_0, k, i],
                    )
                    T = air_temperature(ts)
                    q_liq = PhasePartition(ts).liq
                    q_vap = q_tot - q_liq
                    θ = dry_pottemp(ts)

                    # microphysics rain should apply here to update qt, ql, thetal, T
                    # now collect the inner values of
                    inner_env[i_ql] += q_liq * weights[j_θ] * sqpi_inv
                    inner_env[i_T] += T * weights[j_θ] * sqpi_inv
                    inner_env[i_θl] += θ * weights[j_θ] * sqpi_inv
                    inner_env[i_ρ] += ρ * weights[j_θ] * sqpi_inv
                    # rain area fraction
                    if aux[:qr_src, k, i] > 0
                        inner_env[i_rf] += weights[j_θ] * sqpi_inv
                    end
                    # cloudy/dry categories for buoyancy in TKE
                    if q_liq > 0.0
                        inner_env[i_cf] += weights[j_θ] * sqpi_inv
                        inner_env[i_qt_cld] += qt_hat * weights[j_θ] * sqpi_inv
                        inner_env[i_T_cld] += T * weights[j_θ] * sqpi_inv
                    else
                        inner_env[i_qt_dry] += qt_hat * weights[j_θ] * sqpi_inv
                        inner_env[i_T_dry] += T * weights[j_θ] * sqpi_inv
                    end
                    # products for variance and covariance source terms
                    inner_src[i_Sqt] += -qr_src * weights[j_θ] * sqpi_inv
                    inner_src[i_SH] += θl_rain_src * weights[j_θ] * sqpi_inv
                    inner_src[i_Sqt_H] += -qr_src * θ * weights[j_θ] * sqpi_inv
                    inner_src[i_SH_H] +=
                        θl_rain_src * θ * weights[j_θ] * sqpi_inv
                    inner_src[i_Sqt_qt] +=
                        -qr_src * q_tot * weights[j_θ] * sqpi_inv
                    inner_src[i_SH_qt] +=
                        θl_rain_src * q_tot * weights[j_θ] * sqpi_inv
                end

                outer_env[i_ql] += inner_env[i_ql] * weights[j_qt] * sqpi_inv
                outer_env[i_T] += inner_env[i_T] * weights[j_qt] * sqpi_inv
                outer_env[i_θl] += inner_env[i_θl] * weights[j_qt] * sqpi_inv
                outer_env[i_ρ] += inner_env[i_ρ] * weights[j_qt] * sqpi_inv
                outer_env[i_cf] += inner_env[i_cf] * weights[j_qt] * sqpi_inv
                outer_env[i_qt_cld] +=
                    inner_env[i_qt_cld] * weights[j_qt] * sqpi_inv
                outer_env[i_qt_dry] +=
                    inner_env[i_qt_dry] * weights[j_qt] * sqpi_inv
                outer_env[i_T_cld] +=
                    inner_env[i_T_cld] * weights[j_qt] * sqpi_inv
                outer_env[i_T_dry] +=
                    inner_env[i_T_dry] * weights[j_qt] * sqpi_inv
                outer_env[i_rf] += inner_env[i_rf] * weights[j_qt] * sqpi_inv

                outer_src[i_Sqt] += inner_src[i_Sqt] * weights[j_qt] * sqpi_inv
                outer_src[i_SH] += inner_src[i_SH] * weights[j_qt] * sqpi_inv
                outer_src[i_Sqt_H] +=
                    inner_src[i_Sqt_H] * weights[j_qt] * sqpi_inv
                outer_src[i_SH_H] +=
                    inner_src[i_SH_H] * weights[j_qt] * sqpi_inv
                outer_src[i_Sqt_qt] +=
                    inner_src[i_Sqt_qt] * weights[j_qt] * sqpi_inv
                outer_src[i_SH_qt] +=
                    inner_src[i_SH_qt] * weights[j_qt] * sqpi_inv
            end

            # aux[:CF, k, i]           = outer_env[i_cf]
            # aux[:t_cloudy, k, i]     = outer_env[i_T]
            # aux[:q_tot_cloudy, k, i] = outer_env[i_qt_cld]
            # aux[:q_vap_cloudy, k, i] = outer_env[i_qt_cld] - outer_env[i_ql]
            # aux[:q_tot_dry, k, i]    = outer_env[i_qt_dry]
            # ts = TemperatureSHumEquil(param_set, aux[:t_cloudy, k], p_0=aux[:p_0, k, i], aux[:q_tot_cloudy, k, i])
            # aux[:θ_cloudy, k, i] = liquid_ice_pottemp(ts)
            # aux[:θ_dry, k, i] = dry_pottemp(ts)

            aux[:CF, k] = outer_env[i_cf]
            aux[:t_cloudy, k] = outer_env[i_T]
            aux[:q_tot_cloudy, k] = outer_env[i_qt_cld]
            aux[:q_vap_cloudy, k] = outer_env[i_qt_cld] - outer_env[i_ql]
            aux[:q_tot_dry, k] = outer_env[i_qt_dry]
            ts = TemperatureSHumEquil(
                param_set,
                aux[:t_cloudy, k],
                p_0 = aux[:p_0, k, i],
                aux[:q_tot_cloudy, k, i],
            )
            aux[:θ_cloudy, k] = liquid_ice_pottemp(ts)
            aux[:θ_dry, k] = dry_pottemp(ts)
        end
    end
end

function compute_subdomain_statistics!(
    grid::Grid{FT},
    i,
    q,
    aux,
    params,
    model::LogNormalQuadrature,
) where {FT, IT}
    gm, en, ud, sd, al = allcombinations(q)
    abscissas, weights = gausshermite(model.order)
    sqrt2 = sqrt(2.0)
    sqpi_inv = 1 / sqrt(pi)
    env_len = 10
    src_len = 6
    @unpack param_set = params

    inner_env = zeros(env_len)
    outer_env = zeros(env_len)
    inner_src = zeros(src_len)
    outer_src = zeros(src_len)
    i_Sqt = collect(1:src_len)
    i_SH = collect(1:src_len)
    i_Sqt_H = collect(1:src_len)
    i_Sqt_qt = collect(1:src_len)
    i_SH_H = collect(1:src_len)
    i_SH_qt = collect(1:src_len)
    i_ql = collect(1:env_len)
    i_T = collect(1:env_len)
    i_θl = collect(1:env_len)
    i_ρ = collect(1:env_len)
    i_cf = collect(1:env_len)
    i_qt_cld = collect(1:env_len)
    i_qt_dry = collect(1:env_len)
    i_T_cld = collect(1:env_len)
    i_T_dry = collect(1:env_len)
    i_rf = collect(1:env_len)

    @inbounds for k in over_elems_real(grid)
        if (
            q[:cv_q_tot, k, i] > eps(FT) &&
            q[:cv_θ_liq, k, i] > eps(FT) &&
            fabs(q[:cv_θ_liq_q_tot, k, i]) > eps(FT) &&
            q[:q_tot, k, i] > eps(FT) &&
            sqrt(q[:cv_q_tot, k, i]) < q[:q_tot, k, i]
        )
            σ_q = sqrt(
                log(q[:cv_q_tot, k, i] / q[:q_tot, k, i] / q[:q_tot, k, i]) + 1,
            )
            σ_h = sqrt(
                log(q[:cv_θ_liq, k, i] / q[:θ_liq, k, i] / q[:θ_liq, k, i]) + 1,
            )
            # Enforce Schwarz's inequality
            corr = max(
                min(q[:cv_θ_liq_q_tot, k, i] / max(σ_h * σ_q, 1e-13), 1),
                -1,
            )
            sd2_hq = log(
                corr * sqrt(q[:cv_q_tot, k, i] * q[:cv_θ_liq, k, i]) /
                (q[:θl, k, i] * q[:q_tot, k, i]) + 1.0,
            )
            σ_cond_θ_q = sqrt(max(σ_h * σ_h - sd2_hq * sd2_hq / σ_q / σ_q, 0.0))
            μ_q = log(
                q[:q_tot, k, i] * q[:q_tot, k, i] / sqrt(
                    q[:q_tot, k, i] * q[:q_tot, k, i] + q[:cv_q_tot, k, i],
                ),
            )
            μ_θ = log(
                q[:θ_liq, k, i] * q[:θ_liq, k, i] / sqrt(
                    q[:θ_liq, k, i] * q[:θ_liq, k, i] + q[:cv_θ_liq, k, i],
                ),
            )
            # clean outer vectors
            outer_src = 0.0 * outer_src
            outer_env = 0.0 * outer_env

            for j_qt in 1:(model.order)
                qt_hat = exp(μ_q + sqrt2 * σ_q * abscissas[j_q])
                μ_θ_star = μ_θ + sd2_hq / sd_q / sd_q * (log(qt_hat) - μ_q)
                # clean innner vectors
                inner_src = 0.0 * inner_src
                inner_env = 0.0 * inner_env
                for j_θ in 1:(model.order)
                    θ_hat = exp(μ_θ_star + sqrt2 * σ_cond_θ_q * abscissas[j_θ])
                    ts = DummyThermoState(
                        param_set;
                        θ = θ_hat,
                        q_tot = qt_hat,
                        ρ_0 = aux[:ρ_0, k, i],
                        p_0 = aux[:p_0, k, i],
                    )
                    T = air_temperature(ts)
                    q_liq = PhasePartition(ts).liq
                    q_vap = q_tot - q_liq
                    θ = dry_pottemp(ts)

                    # microphysics rain should apply here to update qt, ql, thetal, T
                    # now collect the inner values of
                    inner_env[i_ql] += q_liq * weights[j_θ] * sqpi_inv
                    inner_env[i_T] += T * weights[j_θ] * sqpi_inv
                    inner_env[i_θl] += θ * weights[j_θ] * sqpi_inv
                    inner_env[i_ρ] += ρ * weights[j_θ] * sqpi_inv
                    # rain area fraction
                    if aux[:qr_src, k, i] > 0
                        inner_env[i_rf] += weights[j_θ] * sqpi_inv
                    end
                    # cloudy/dry categories for buoyancy in TKE
                    if q_liq > 0.0
                        inner_env[i_cf] += weights[j_θ] * sqpi_inv
                        inner_env[i_qt_cld] += qt_hat * weights[j_θ] * sqpi_inv
                        inner_env[i_T_cld] += T * weights[j_θ] * sqpi_inv
                    else
                        inner_env[i_qt_dry] += qt_hat * weights[j_θ] * sqpi_inv
                        inner_env[i_T_dry] += T * weights[j_θ] * sqpi_inv
                    end
                    # products for variance and covariance source terms
                    inner_src[i_Sqt] += -qr_src * weights[j_θ] * sqpi_inv
                    inner_src[i_SH] += θl_rain_src * weights[j_θ] * sqpi_inv
                    inner_src[i_Sqt_H] += -qr_src * θ * weights[j_θ] * sqpi_inv
                    inner_src[i_SH_H] +=
                        θl_rain_src * θ * weights[j_θ] * sqpi_inv
                    inner_src[i_Sqt_qt] +=
                        -qr_src * q_tot * weights[j_θ] * sqpi_inv
                    inner_src[i_SH_qt] +=
                        θl_rain_src * q_tot * weights[j_θ] * sqpi_inv
                end

                outer_env[i_ql] += inner_env[i_ql] * weights[j_qt] * sqpi_inv
                outer_env[i_T] += inner_env[i_T] * weights[j_qt] * sqpi_inv
                outer_env[i_θl] += inner_env[i_θl] * weights[j_qt] * sqpi_inv
                outer_env[i_ρ] += inner_env[i_ρ] * weights[j_qt] * sqpi_inv
                outer_env[i_cf] += inner_env[i_cf] * weights[j_qt] * sqpi_inv
                outer_env[i_qt_cld] +=
                    inner_env[i_qt_cld] * weights[j_qt] * sqpi_inv
                outer_env[i_qt_dry] +=
                    inner_env[i_qt_dry] * weights[j_qt] * sqpi_inv
                outer_env[i_T_cld] +=
                    inner_env[i_T_cld] * weights[j_qt] * sqpi_inv
                outer_env[i_T_dry] +=
                    inner_env[i_T_dry] * weights[j_qt] * sqpi_inv
                outer_env[i_rf] += inner_env[i_rf] * weights[j_qt] * sqpi_inv

                outer_src[i_Sqt] += inner_src[i_Sqt] * weights[j_qt] * sqpi_inv
                outer_src[i_SH] += inner_src[i_SH] * weights[j_qt] * sqpi_inv
                outer_src[i_Sqt_H] +=
                    inner_src[i_Sqt_H] * weights[j_qt] * sqpi_inv
                outer_src[i_SH_H] +=
                    inner_src[i_SH_H] * weights[j_qt] * sqpi_inv
                outer_src[i_Sqt_qt] +=
                    inner_src[i_Sqt_qt] * weights[j_qt] * sqpi_inv
                outer_src[i_SH_qt] +=
                    inner_src[i_SH_qt] * weights[j_qt] * sqpi_inv
            end

            aux[:CF, k, i] = outer_env[i_cf]
            aux[:t_cloudy, k, i] = outer_env[i_T]
            aux[:q_tot_cloudy, k, i] = outer_env[i_qt_cld]
            aux[:q_vap_cloudy, k, i] = outer_env[i_qt_cld] - outer_env[i_ql]
            aux[:q_tot_dry, k, i] = outer_env[i_qt_dry]
            ts = TemperatureSHumEquil(
                param_set,
                aux[:t_cloudy, k],
                p_0 = aux[:p_0, k, i],
                aux[:q_tot_cloudy, k, i],
            )
            aux[:θ_cloudy, k, i] = liquid_ice_pottemp(ts)
            aux[:θ_dry, k, i] = dry_pottemp(ts)

        end
    end
end
