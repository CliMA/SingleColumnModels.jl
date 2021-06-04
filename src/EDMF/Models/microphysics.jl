##### Microphysics models

"""
    compute_cloud_phys!

Define microphysics fields
 - `tmp[:CF, k, i]`
 - `tmp[:θ_cloudy, k, i]`
 - `tmp[:t_cloudy, k, i]`
 - `tmp[:q_tot_cloudy, k, i]`
 - `tmp[:q_vap_cloudy, k, i]`
 - `tmp[:θ_dry, k, i]`
 - `tmp[:q_tot_dry, k, i]`

 for all `k` and all `i`
"""
function compute_cloud_phys! end

function compute_cloud_phys!(grid::Grid{FT}, q, tmp, params) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @unpack params param_set
    @inbounds for k in over_elems_real(grid)
        q_tot = q[:q_tot, k, en]
        ts = ActiveThermoState(param_set, q, tmp, k, en)
        T = air_temperature(ts)
        q_liq = PhasePartition(ts).liq
        q_vap = q_tot - q_liq
        θ = dry_pottemp(ts)
        if q_liq > 0
            tmp[:CF, k] = FT(1)
            tmp[:θ_cloudy, k] = θ
            tmp[:t_cloudy, k] = T
            tmp[:q_tot_cloudy, k] = q_tot
            tmp[:q_vap_cloudy, k] = q_vap
        else
            tmp[:CF, k] = FT(0)
            tmp[:θ_dry, k] = θ
            tmp[:q_tot_dry, k] = q_tot
        end
    end
end
