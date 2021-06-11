function compute_new_ud_a!(grid, q_new, q, q_tendencies, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            a_predict = q[:a, k, i] + params[:Δt][1] * q_tendencies[:a, k, i]
            q_new[:a, k, i] = bound(a_predict, params[:a_bounds])
        end
    end
end

function compute_new_ud_w!(grid, q_new, q, q_tendencies, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    # Solve for updraft velocity
    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            a_new_k = q_new[:a, k, i]
            ρ_k = aux[:ρ_0, k]
            w_i = q[:w, k, i]
            a_k = q[:a, k, i]
            ρa_k = ρ_k * a_k
            ρa_new_k = ρ_k * a_new_k
            ρaw_k = ρa_k * w_i
            w_predict =
                ρaw_k / ρa_new_k +
                params[:Δt][1] / ρa_new_k * q_tendencies[:w, k, i]
            q_new[:w, k, i] = bound(w_predict, params[:w_bounds])
        end
    end
end

function compute_new_ud_scalars!(grid, q_new, q, q_tendencies, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            a_k = q[:a, k, i]
            a_k_new = q_new[:a, k, i]
            ρ_k = aux[:ρ_0, k]
            ρa_k = ρ_k * a_k
            ρa_new_k = ρ_k * a_k_new
            θ_liq_predict =
                (
                    ρa_k * q[:θ_liq, k, i] +
                    params[:Δt][1] * q_tendencies[:θ_liq, k, i]
                ) / ρa_new_k
            q_tot_predict =
                (
                    ρa_k * q[:q_tot, k, i] +
                    params[:Δt][1] * q_tendencies[:q_tot, k, i]
                ) / ρa_new_k
            q_new[:θ_liq, k, i] = θ_liq_predict
            q_new[:q_tot, k, i] = q_tot_predict
        end
    end
end

#####
##### Update functions
#####

function assign_new_to_values!(grid, q_new, q, aux)
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for k in over_elems(grid), i in ud, name in (:w, :q_tot, :θ_liq)
        q_new[name, k, i] = q[name, k, i]
    end
end

function assign_values_to_new!(grid, q, q_new, aux)
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for k in over_elems(grid), i in ud
        @inbounds for name in (:a, :w, :q_tot, :θ_liq)
            q[name, k, i] = q_new[name, k, i]
        end
    end
    @inbounds for k in over_elems(grid)
        q[:tke, k, en] = q_new[:tke, k, en]
        q[:a, k, en] = 1 - sum([q_new[:a, k, i] for i in ud])
        q[:θ_liq, k, gm] = q_new[:θ_liq, k, gm]
        q[:q_tot, k, gm] = q_new[:q_tot, k, gm]
    end
end
