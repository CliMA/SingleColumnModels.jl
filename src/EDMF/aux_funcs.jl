#### AuxiliaryFuncs

export ActiveThermoState

export update_dt!
export heaviside

"""
    heaviside(x_1, x_2)

Heaviside function
 - `x_1`
 - `x_2` value specified at edge case `x_1 = 0`
"""
heaviside(x_1, x_2) = x_1 == 0 ? x_2 : typeof(x_1)(x_1 > 0)

"""
    ActiveThermoState(q, tmp, k, i)

Returns a `ThermodynamicState` using grid-mean
quantities at element `k`.
"""
@inline function ActiveThermoState(param_set, q, tmp, k, i)
    return LiquidIcePotTempSHumEquil_old(
        param_set,
        q[:θ_liq, k, i],
        q[:q_tot, k, i],
        tmp[:ρ_0, k],
        tmp[:p_0, k],
    )
end

"""
    ActiveThermoState(q, tmp, k, i)

Returns a `ThermodynamicState` using grid-mean
quantities at element `k`.
"""
@inline function ActiveThermoState(param_set, q, tmp, k::Dual, i)
    return LiquidIcePotTempSHumEquil_old.(
        Ref(param_set),
        q[:θ_liq, k, i],
        q[:q_tot, k, i],
        tmp[:ρ_0, k],
        tmp[:p_0, k],
    )
end

function update_dt!(grid, params, q, t)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack params Δt Δt_min
    u_max = max([q[:w, k, i] for i in ud for k in over_elems(grid)]...)
    Δt = [min(Δt_min, 0.5 * grid.Δz / max(u_max, 1e-10))]
    Δti = [1 / Δt[1]]
    params[:Δt] = Δt
    params[:Δti] = Δti

    # println("----------------------------")
    t[1] += Δt[1]
    percent_done = t[1] / params[:t_end] * 100.0
    # @show t[1], Δt[1], percent_done, params[:t_end]
end
