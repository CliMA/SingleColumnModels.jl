#### SolveTurbConv

function solve!(
    q_new::StateVec,
    q::StateVec,
    q_tendencies::StateVec,
    ∑tendencies_params,
)
    @unpack params, grid, aux_O2, aux, tri_diag, case = ∑tendencies_params
    gm, en, ud, sd, al = allcombinations(q)

    compute_new_ud_a!(grid, q_new, q, q_tendencies, aux, params)
    apply_bcs!(grid, q_new, aux, params, case)

    compute_new_ud_w!(grid, q_new, q, q_tendencies, aux, params)
    compute_new_ud_scalars!(grid, q_new, q, q_tendencies, aux, params)

    apply_bcs!(grid, q_new, aux, params, case)

    solve_tridiag_wrapper!(grid, q_new, :tke, en, tri_diag[:tke])
    solve_tridiag_wrapper!(grid, q_new, :q_tot, gm, tri_diag[:q_tot])
    solve_tridiag_wrapper!(grid, q_new, :θ_liq, gm, tri_diag[:θ_liq])

    assign_values_to_new!(grid, q, q_new, aux)
    apply_bcs!(grid, q, aux, params, case)
    extrap_0th_order!(q, (:θ_liq, :q_tot), grid, gm)

end
