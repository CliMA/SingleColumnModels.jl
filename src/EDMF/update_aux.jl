#### PrecomputeVars

function pre_compute_vars!(grid, q, aux, aux_O2, UpdVar, params)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set = params

    diagnose_environment!(q, grid, :a, (:q_tot, :θ_liq, :w))

    saturation_adjustment_sd!(grid, q, aux, params)

    @inbounds for k in over_elems_real(grid)
        ts = ActiveThermoState(param_set, q, aux, k, gm)
        aux[:θ_ρ, k] = virtual_pottemp(ts)
    end
    params[:zi] = compute_inversion_height(aux, q, grid, params)
    params[:wstar] = compute_convective_velocity(params[:bflux], params[:zi])

    compute_entrainment_detrainment!(
        grid,
        UpdVar,
        aux,
        q,
        params,
        params[:EntrDetrModel],
    )
    compute_cloud_phys!(grid, q, aux, params)
    compute_buoyancy!(grid, q, aux, params)
    compute_pressure!(grid, q, aux, params, params[:PressureModel])

    filter_scalars!(grid, q, aux, params)

    compute_cv_gm!(grid, q, :w, :w, :tke, 0.5)
    compute_mf_gm!(grid, q, aux)
    compute_mixing_length!(grid, q, aux, params, params[:MixingLengthModel])
    compute_eddy_diffusivities_tke!(
        grid,
        q,
        aux,
        params,
        params[:EddyDiffusivityModel],
    )

    compute_tke_buoy!(grid, q, aux, aux_O2, :tke, params)
    compute_cv_entr!(grid, q, aux, aux_O2, :w, :w, :tke, 0.5)
    compute_cv_shear!(grid, q, aux, aux_O2, :w, :w, :tke)
    compute_cv_interdomain_src!(grid, q, aux, aux_O2, :w, :w, :tke, 0.5)
    compute_tke_pressure!(
        grid,
        q,
        aux,
        aux_O2,
        :tke,
        params,
        params[:PressureModel],
    )
    compute_cv_env!(grid, q, aux, aux_O2, :w, :w, :tke, 0.5)

    cleanup_covariance!(grid, q)

end
