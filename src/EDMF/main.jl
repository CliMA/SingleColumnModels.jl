#### Main

"""
    run(case)
Solve the Eddy-Diffusivity Mass-Flux (EDMF) equations for a
stand-alone `case`
"""
function run(param_set, case)
    params = Params(param_set, case)

    tc = TurbConv(params, case)

    grid = tc.grid
    q = tc.q
    q_new = tc.q_new
    tmp = tc.tmp
    q_tendencies = tc.tendencies
    dir_tree = tc.dir_tree
    tri_diag = tc.tri_diag
    tmp_O2 = tc.tmp_O2

    gm, en, ud, sd, al = allcombinations(tc.q)
    init_ref_state!(tmp, grid, params, dir_tree)
    init_forcing!(q, tmp, grid, params, dir_tree, case)
    init_state_vecs!(q, tmp, grid, params, dir_tree, case)


    params[:UpdVar] =
        [UpdraftVar(0, params[:SurfaceModel].area, length(ud)) for i in al]
    # export_initial_conditions(q, tmp, grid, dir_tree[:processed_initial_conditions], true)

    @unpack params Δt t_end

    i_Δt, i_export, t = [0], [0], [0.0]

    assign!(q_tendencies, (:u, :v, :q_tot, :θ_liq), grid, 0.0)
    update_surface!(tmp, q, grid, params, params[:SurfaceModel])
    update_forcing!(q_tendencies, tmp, q, grid, params, params[:ForcingType])
    compute_cloud_base_top_cover!(params[:UpdVar], grid, q, tmp)

    pre_compute_vars!(grid, q, tmp, tmp_O2, params[:UpdVar], params)

    apply_bcs!(grid, q, tmp, params, case)

    nc_q = NetCDFWriter(joinpath(dir_tree[:solution_raw], "prog_vs_time"))
    nc_q = init_data(nc_q, grid, q)
    nc_tmp = NetCDFWriter(joinpath(dir_tree[:solution_raw], "aux_vs_time"))
    nc_tmp = init_data(nc_tmp, grid, tmp)

    @show joinpath(pwd(), full_name(nc_q))

    while t[1] < t_end
        assign!(q_tendencies, (:u, :v, :q_tot, :θ_liq), grid, 0.0)

        update_surface!(tmp, q, grid, params, params[:SurfaceModel])
        update_forcing!(
            q_tendencies,
            tmp,
            q,
            grid,
            params,
            params[:ForcingType],
        )

        pre_compute_vars!(grid, q, tmp, tmp_O2, params[:UpdVar], params)

        update!(
            grid,
            q_new,
            q,
            q_tendencies,
            tmp,
            tmp_O2,
            case,
            tri_diag,
            params,
        )

        update_dt!(grid, params, q, t)

        grid_mean!(tmp, q, :a, (:T, :q_liq, :buoy), grid)
        compute_cloud_base_top_cover!(params[:UpdVar], grid, q, tmp)

        i_export[1] += 1
        if mod(i_export[1], params[:export_frequency]) == 0
            println("********* EXPORTING *********")
            @show nc_q.filename
            append_data(nc_q, q, t[1])
            append_data(nc_tmp, tmp, t[1])
        end
    end
    extrap_0th_order!(q, (:θ_liq, :q_tot), grid, gm)
    extrap_0th_order!(tmp, :T, grid, gm)

    append_data(nc_q, q, t[1])
    append_data(nc_tmp, tmp, t[1])

    return (grid, q, tmp, params)
end
