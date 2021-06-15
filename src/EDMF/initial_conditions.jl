#### InitialConditions

const Rd_SCAMPY = 287.1
const eps_vi = 1.60745384883
const dycoms_cp = 1.60745384883

"""
    init_state_vecs!

Defines initial conditions for state vectors `q` and `aux` for all sub-domains.
"""
function init_state_vecs! end

"""
    initialize_updrafts!

Defines initial conditions for state vectors `q` and `aux` for the updraft sub-domains.
"""
function initialize_updrafts! end

"""
    init_forcing!

Defines initial conditions for forcing terms in the state vector `aux` for all sub-domains.
"""
function init_forcing! end


function initialize_updrafts!(
    q::StateVec,
    aux::StateVec,
    grid::Grid,
    params,
    output_dir::String,
    ::BOMEX,
)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack SurfaceModel = params
    k_1 = first_interior(grid, Zmin())
    n_updrafts = length(ud)
    for i in ud
        for k in over_elems(grid)
            q[:w, k, i] = 0
            q[:a, k, i] = bound(0.0, params[:a_bounds])
        end
        q[:a, k_1, i] = bound(SurfaceModel.area / n_updrafts, params[:a_bounds])
    end
    return
end

function init_state_vecs!(
    q::StateVec,
    aux::StateVec,
    grid::Grid,
    params,
    output_dir::String,
    case::BOMEX,
)
    @unpack a_bounds, SurfaceModel, param_set = params
    z = grid.zc

    gm, en, ud, sd, al = allcombinations(q)
    k_1 = first_interior(grid, Zmin())

    @inbounds for k in over_elems(grid)
        @inbounds for ϕ in var_names(q)
            @inbounds for i in over_sub_domains(q, ϕ)
                q[ϕ, k, i] = 0.0
            end
        end

        q[:a, k, gm] = 1.0
        for i in ud
            q[:a, k, i] = bound(0.0, a_bounds)
            k == k_1 && (q[:a, k, i] = bound(SurfaceModel.area, a_bounds))
        end
        q[:a, k, en] = q[:a, k, gm] - sum([q[:a, k, i] for i in ud])

        # Set qt profile
        if z[k] <= 520
            q[:q_tot, k, gm] = (17.0 + (z[k]) * (16.3 - 17.0) / 520.0) / 1000.0
        end
        if z[k] > 520.0 && z[k] <= 1480.0
            q[:q_tot, k, gm] =
                (16.3 + (z[k] - 520.0) * (10.7 - 16.3) / (1480.0 - 520.0)) /
                1000.0
        end
        if z[k] > 1480.0 && z[k] <= 2000.0
            q[:q_tot, k, gm] =
                (10.7 + (z[k] - 1480.0) * (4.2 - 10.7) / (2000.0 - 1480.0)) /
                1000.0
        end
        if z[k] > 2000.0
            q[:q_tot, k, gm] =
                (4.2 + (z[k] - 2000.0) * (3.0 - 4.2) / (3000.0 - 2000.0)) /
                1000.0
        end

        # Set θ_liq profile
        if z[k] <= 520.0
            q[:θ_liq, k, gm] = 298.7
        end
        if z[k] > 520.0 && z[k] <= 1480.0
            q[:θ_liq, k, gm] =
                298.7 + (z[k] - 520) * (302.4 - 298.7) / (1480.0 - 520.0)
        end
        if z[k] > 1480.0 && z[k] <= 2000
            q[:θ_liq, k, gm] =
                302.4 + (z[k] - 1480.0) * (308.2 - 302.4) / (2000.0 - 1480.0)
        end
        if z[k] > 2000.0
            q[:θ_liq, k, gm] =
                308.2 + (z[k] - 2000.0) * (311.85 - 308.2) / (3000.0 - 2000.0)
        end
    end

    @inbounds for k in over_elems_real(grid)
        ts = ActiveThermoState(param_set, q, aux, k, gm)

        aux[:T, k, gm] = air_temperature(ts)

        # Set u profile
        if z[k] <= 700.0
            q[:u, k, gm] = -8.75
        end
        if z[k] > 700.0
            q[:u, k, gm] =
                -8.75 + (z[k] - 700.0) * (-4.61 - -8.75) / (3000.0 - 700.0)
        end
    end # end over_elems

    # Extrapolate to ghost points
    extrap_0th_order!(q, (:θ_liq, :q_tot, :u), grid, gm)
    extrap_0th_order!(aux, :T, grid, gm)
    # Use grid-mean for sub-domain values:

    initialize_updrafts!(q, aux, grid, params, output_dir, case)
    distribute!(q, grid, (:q_tot, :θ_liq))
    distribute!(aux, grid, (:q_liq, :T))
    diagnose_environment!(q, grid, :a, (:q_tot, :θ_liq, :w))
end

function init_forcing!(
    q::StateVec,
    aux::StateVec,
    grid::Grid{FT},
    params,
    output_dir::String,
    case::BOMEX,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set = params
    z = grid.zc
    for k in over_elems_real(grid)
        # Geostrophic velocity profiles. vg = 0
        aux[:ug, k, gm] = -10.0 + (1.8e-3) * z[k]
        Π = TD.exner_given_pressure(param_set, aux[:p_0, k])
        # Set large-scale cooling
        if z[k] <= 1500.0
            aux[:dTdt, k, gm] = (-2.0 / (3600 * 24.0)) * Π
        else
            aux[:dTdt, k, gm] =
                (
                    -2.0 / (3600 * 24.0) +
                    (z[k] - 1500.0) * (0.0 - -2.0 / (3600 * 24.0)) /
                    (3000.0 - 1500.0)
                ) * Π
        end
        # Set large-scale drying
        if z[k] <= 300.0
            aux[:dqtdt, k, gm] = -1.2e-8   #kg/(kg * s)
        end
        if z[k] > 300.0 && z[k] <= 500.0
            aux[:dqtdt, k, gm] =
                -1.2e-8 + (z[k] - 300.0) * (0.0 - -1.2e-8) / (500.0 - 300.0) #kg/(kg * s)
        end

        #Set large scale subsidence
        if z[k] <= 1500.0
            aux[:subsidence, k, gm] =
                0.0 + z[k] * (-0.65 / 100.0 - 0.0) / (1500.0 - 0.0)
        end
        if z[k] > 1500.0 && z[k] <= 2100.0
            aux[:subsidence, k, gm] =
                -0.65 / 100 +
                (z[k] - 1500.0) * (0.0 - -0.65 / 100.0) / (2100.0 - 1500.0)
        end
    end
end

function init_forcing!(
    q::StateVec,
    aux::StateVec,
    grid::Grid{FT},
    params,
    output_dir::String,
    case::DYCOMS,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set, divergence, F0, F1, kappa, alpha_z = params
    k_1 = first_interior(grid, Zmin())
    z = grid.zc
    for k in over_elems(grid)
        # Geostrophic velocity profiles. vg = 0
        aux[:ug, k, gm] = 7.0
        aux[:vg, k, gm] = -5.5
        # Set large scale subsidence
        aux[:subsidence, k, gm] = -z[k] * divergence
        aux[:dqtdt, k, gm] = 0
    end

    ze = grid.ze
    # see eq. 3 in Stevens et. al. 2005 DYCOMS paper
    # find zi (level of 8.0 g/kg isoline of qt)
    local zi, ρ_i
    for k in over_elems(grid)
        if q[:q_tot, k, gm] < 8.0 / 1000
            idx_zi = k
            # will be used at cell edges
            zi = ze[idx_zi]
            ρ_i = aux[:ρ_0, idx_zi]
            break
        end
    end

    # cloud-top cooling
    q_0 = 0.0

    aux[:f_rad, k_1] = F0 * exp(-q_0)
    for k in reverse(over_elems_real(grid))
        q_0 += kappa * aux[:ρ_0, k] * aux[:q_liq, k] * grid.Δz
        aux[:f_rad, k] = F0 * exp(-q_0)
    end

    # cloud-base warming
    q_1 = 0.0
    aux[:f_rad, k_1] += F1 * exp(-q_1)
    for k in over_elems_real(grid)
        q_1 += kappa * aux[:ρ_0, k - 1] * aux[:q_liq, k - 1] * grid.Δz
        aux[:f_rad, k] += F1 * exp(-q_1)
    end

    # cooling in free troposphere
    for k in over_elems_real(grid)
        if z[k] > zi
            cbrt_z = cbrt(z[k] - zi)
            aux[:f_rad, k] +=
                ρ_i * dycoms_cp * divergence * alpha_z * (cbrt_z^4) / 4 +
                zi * cbrt_z
        end
    end
    # condition at the top
    k_n = first_interior(grid, Zmax())
    cbrt_z = cbrt(z[k_n] + grid.Δz - zi)
    aux[:f_rad, k_n] +=
        ρ_i * dycoms_cp * divergence * alpha_z * (cbrt_z^4 / 4 + zi * cbrt_z)
    extrap_0th_order!(aux, (:f_rad,), grid, gm)

    for k in over_elems_real(grid)
        aux[:dTdt, k] =
            -(aux[:f_rad, k + 1] - aux[:f_rad, k]) / grid.Δz / aux[:ρ_0, k] /
            dycoms_cp
    end


    nc = NetCDFWriter(joinpath(output_dir, "ic_forcing"))
    export_state(nc, grid, aux)
end

rho_c(p0, T, qt, qv) = p0 / ((Rd_SCAMPY * T) * (1.0 - qt + eps_vi * qv))

buoyancy_c(param_set, rho0::FT, rho::FT) where {FT} =
    grav(param_set) * (rho0 - rho) / rho0

function init_state_vecs!(
    q::StateVec,
    aux::StateVec,
    grid::Grid,
    params,
    output_dir::String,
    case::DYCOMS,
)
    @unpack a_bounds, SurfaceModel, param_set = params
    z = grid.zc

    gm, en, ud, sd, al = allcombinations(q)

    @inbounds for k in over_elems(grid)
        @inbounds for ϕ in var_names(q)
            @inbounds for i in over_sub_domains(q, ϕ)
                if z[k] <= 840.0
                    q[:θ_liq, k] = 289.0
                end
                if z[k] > 840.0
                    q[:θ_liq, k] = (297.5 + (z[k] - 840.0)^(1.0 / 3.0))
                end

                # qt profile as defined in DYCOMS
                if z[k] <= 840.0
                    q[:q_tot, k] = 9.0 / 1000.0
                end
                if z[k] > 840.0
                    q[:q_tot, k] = 1.5 / 1000.0
                end

                # ql and T profile
                # (calculated by saturation adjustment using thetal and qt values provided in DYCOMS
                # and using Rd, cp and L constants as defined in DYCOMS)
                # aux[:T, k], aux[:q_liq, k] = dycoms_sat_adjst(aux[:p_0, k], q[:θ_liq, k], q[:q_tot, k])
                p_0 = aux[:p_0, k]
                ρ_0 = aux[:ρ_0, k]
                θ_liq_ice = q[:θ_liq, k]
                q_tot = q[:q_tot, k]
                ts = LiquidIcePotTempSHumEquil_old(
                    param_set,
                    θ_liq_ice,
                    q_tot,
                    ρ_0,
                    p_0,
                )
                aux[:T, k] = air_temperature(ts)
                aux[:q_liq, k] = liquid_specific_humidity(ts)

                # buoyancy profile
                qv = q[:q_tot, k] - aux[:q_ice, k] - aux[:q_liq, k]
                rho = rho_c(aux[:p_0, k], aux[:T, k], q[:q_tot, k], qv)
                aux[:buoy, k] = buoyancy_c(param_set, aux[:p_0, k], rho)

                # velocity profile (geostrophic)
                aux[:ug, k] = 7.0
                aux[:vg, k] = -5.5

            end
        end
    end

end
