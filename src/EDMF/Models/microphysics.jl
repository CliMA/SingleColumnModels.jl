##### Microphysics models

"""
    compute_cloud_phys!

Define microphysics fields
 - `aux[:CF, k, i]`
 - `aux[:θ_cloudy, k, i]`
 - `aux[:t_cloudy, k, i]`
 - `aux[:q_tot_cloudy, k, i]`
 - `aux[:q_vap_cloudy, k, i]`
 - `aux[:θ_dry, k, i]`
 - `aux[:q_tot_dry, k, i]`

 for all `k` and all `i`
"""
function compute_cloud_phys! end

function compute_cloud_phys!(grid::Grid{FT}, q, aux, params) where {FT}
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
