#### Update forcing

####
#### Types
####

abstract type ForcingType end

struct NoForcing <: ForcingType end

struct StandardForcing{FT} <: ForcingType
    apply_subsidence::Bool
    apply_coriolis::Bool
    coriolis_param::FT
    StandardForcing(;
        apply_subsidence,
        apply_coriolis,
        coriolis_param::FT,
    ) where {FT} = new{FT}(apply_subsidence, apply_coriolis, coriolis_param)
end

####
#### Helper functions
####

function coriolis_force!(aux, grid, q, coriolis_param)
    gm, en, ud, sd, al = allcombinations(q)
    for k in over_elems_real(grid)
        aux[:forcing_u, k, gm] -= coriolis_param * (aux[:vg, k] - q[:v, k, gm])
        aux[:forcing_v, k, gm] += coriolis_param * (aux[:ug, k] - q[:u, k, gm])
    end
end

####
#### Update forcing methods
####

"""
    update_forcing!

Define forcing fields in the `aux` state.
"""
function update_forcing! end

function update_forcing!(
    aux::StateVec,
    q::StateVec,
    grid::Grid,
    params,
    ::NoForcing,
) end

function update_forcing!(
    aux::StateVec,
    q::StateVec,
    grid::Grid,
    params,
    forcing::StandardForcing,
)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set = params

    for k in over_elems_real(grid)
        # Apply large-scale horizontal advection tendencies
        aux[:forcing_θ_liq, k, gm] +=
            aux[:dTdt, k] / TD.exner_given_pressure(param_set, aux[:p_0, k])
        aux[:forcing_q_tot, k, gm] += aux[:dqtdt, k]
    end
    if forcing.apply_subsidence
        for k in over_elems_real(grid)
            # Apply large-scale subsidence tendencies
            aux[:forcing_θ_liq, k, gm] -=
                grad(q[:θ_liq, Dual(k), gm], grid) * aux[:subsidence, k]
            aux[:forcing_q_tot, k, gm] -=
                grad(q[:q_tot, Dual(k), gm], grid) * aux[:subsidence, k]
        end
    end

    if forcing.apply_coriolis
        coriolis_force!(aux, grid, q, forcing.coriolis_param)
    end
end
