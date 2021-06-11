#### AuxiliaryFuncs

export ActiveThermoState

using LambertW
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
    ActiveThermoState(q, aux, k, i)

Returns a `ThermodynamicState` using grid-mean
quantities at element `k`.
"""
@inline function ActiveThermoState(param_set, q, aux, k, i)
    return PhaseEquil_pθq(
        param_set,
        aux[:p_0, k],
        q[:θ_liq, k, i],
        q[:q_tot, k, i],
    )
end

"""
    ActiveThermoState(q, aux, k, i)

Returns a `ThermodynamicState` using grid-mean
quantities at element `k`.
"""
@inline function ActiveThermoState(param_set, q, aux, k::Dual, i)
    return PhaseEquil_pθq.(
        Ref(param_set),
        aux[:p_0, k],
        q[:θ_liq, k, i],
        q[:q_tot, k, i],
    )
end

function update_dt!(grid, params, q, t)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack Δt, Δt_min = params
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


function lamb_smooth_minimum(l, lower_bound, upper_bound)
    leng = size(l)
    xmin = minimum(l)
    lambda0 = max(
        xmin * lower_bound / real(LambertW.lambertw(2.0 / MathConstants.e)),
        upper_bound,
    )

    i = 1
    num = 0
    den = 0
    while (tuple(i) < leng)
        num += l[i] * exp(-(l[i] - xmin) / lambda0)
        den += exp(-(l[i] - xmin) / lambda0)
        i += 1
    end
    smin = num / den

    return smin
end
