function ∑tendencies!(q_tendencies, q, ∑tendencies_params, t)
    @unpack params, grid, aux_O2, aux, case = ∑tendencies_params

    gm, en, ud, sd, al = allcombinations(q)

    for k in over_elems(grid)
        for (tend, force) in (
            (:u, :u_forcing),
            (:v, :v_forcing),
            (:θ_liq, :θ_liq_forcing),
            (:q_tot, :q_tot_forcing),
        )
            q_tendencies[tend, k] = aux[force, k]
        end
    end

    compute_tendencies_en_O2!(grid, q_tendencies, aux_O2, :tke)
    compute_tendencies_gm_scalars!(grid, q_tendencies, q, aux, params)
    compute_tendencies_ud!(grid, q_tendencies, q, aux, params)
end
