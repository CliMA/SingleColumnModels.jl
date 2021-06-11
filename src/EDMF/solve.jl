#### SolveTurbConv

function solve!(
    q_new::StateVec,
    q::StateVec,
    q_tendencies::StateVec,
    tri_diag::StateVec,
    ∑tendencies_params,
)
    @unpack params, grid, aux_O2, aux, case = ∑tendencies_params
    gm, en, ud, sd, al = allcombinations(q)

    compute_new_ud_a!(grid, q_new, q, q_tendencies, aux, params)
    apply_bcs!(grid, q_new, aux, params, case)

    compute_new_ud_w!(grid, q_new, q, q_tendencies, aux, params)
    compute_new_ud_scalars!(grid, q_new, q, q_tendencies, aux, params)

    apply_bcs!(grid, q_new, aux, params, case)

    compute_new_en_O2!(
        grid,
        q_new,
        q,
        q_tendencies,
        aux,
        aux_O2,
        params,
        :tke,
        tri_diag,
    )

    compute_new_gm_scalars!(grid, q_new, q, q_tendencies, params, aux, tri_diag)

    assign_values_to_new!(grid, q, q_new, aux)
    apply_bcs!(grid, q, aux, params, case)
    extrap_0th_order!(q, (:θ_liq, :q_tot), grid, gm)

end
