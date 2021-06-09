#### SolveTurbConv

function update!(
    grid::Grid,
    q_new::StateVec,
    q::StateVec,
    q_tendencies::StateVec,
    aux::StateVec,
    aux_O2::Dict,
    case::Case,
    tri_diag::StateVec,
    params,
)

    assign_new_to_values!(grid, q_new, q, aux)

    compute_tendencies_en_O2!(grid, q_tendencies, aux_O2, :tke)
    compute_tendencies_gm_scalars!(grid, q_tendencies, q, aux, params)
    compute_tendencies_ud!(grid, q_tendencies, q, aux, params)

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

end
