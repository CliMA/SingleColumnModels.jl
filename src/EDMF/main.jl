#### Main

"""
    run(case)
Solve the Eddy-Diffusivity Mass-Flux (EDMF) equations for a
stand-alone `case`
"""
function run(param_set, case, output_dir)
    params = Params(param_set, case)

    tc = TurbConv(params, case, output_dir)

    grid = tc.grid
    q = tc.q
    q_new = tc.q_new
    aux = tc.aux
    q_tendencies = tc.tendencies
    tri_diag = tc.tri_diag
    aux_O2 = tc.aux_O2

    gm, en, ud, sd, al = allcombinations(tc.q)

    # Initialize nc files
    nc_ics_q = NetCDFWriter(joinpath(output_dir, "ics_q"))
    nc_ics_q = init_data(nc_ics_q, grid, q)
    nc_ics_aux = NetCDFWriter(joinpath(output_dir, "ics_aux"))
    nc_ics_aux = init_data(nc_ics_aux, grid, aux)

    # Initialize reference state
    init_ref_state!(aux, grid, params, output_dir)
    append_data(nc_ics_aux, aux, 0)

    # Initialize forcing
    init_forcing!(q, aux, grid, params, output_dir, case)
    append_data(nc_ics_aux, aux, 0)

    # Initialize solution
    init_state_vecs!(q, aux, grid, params, output_dir, case)
    append_data(nc_ics_q, q, 0)

    params[:UpdVar] =
        [UpdraftVar(0, params[:SurfaceModel].area, length(ud)) for i in al]

    @unpack Δt, t_end = params

    i_Δt, i_export, t = [0], [0], [0.0]

    assign!(q_tendencies, (:u, :v, :q_tot, :θ_liq), grid, 0.0)
    update_surface!(aux, q, grid, params, params[:SurfaceModel])
    update_forcing!(q_tendencies, aux, q, grid, params, params[:ForcingType])
    compute_cloud_base_top_cover!(params[:UpdVar], grid, q, aux)

    update_aux!(grid, q, aux, aux_O2, params[:UpdVar], params)
    # append_data(nc_ics_aux, aux, 0) # will this overwrite anything?

    apply_bcs!(grid, q, aux, params, case)

    nc_q = NetCDFWriter(joinpath(output_dir, "prog_vs_time"))
    nc_q = init_data(nc_q, grid, q)
    nc_aux = NetCDFWriter(joinpath(output_dir, "aux_vs_time"))
    nc_aux = init_data(nc_aux, grid, aux)

    @show joinpath(pwd(), full_name(nc_q))

    while t[1] < t_end
        assign!(q_tendencies, (:u, :v, :q_tot, :θ_liq), grid, 0.0)

        update_surface!(aux, q, grid, params, params[:SurfaceModel])
        update_forcing!(
            q_tendencies,
            aux,
            q,
            grid,
            params,
            params[:ForcingType],
        )

        update_aux!(grid, q, aux, aux_O2, params[:UpdVar], params)

        update!(
            grid,
            q_new,
            q,
            q_tendencies,
            aux,
            aux_O2,
            case,
            tri_diag,
            params,
        )

        update_dt!(grid, params, q, t)

        grid_mean!(aux, q, :a, (:T, :q_liq, :buoy), grid)
        compute_cloud_base_top_cover!(params[:UpdVar], grid, q, aux)

        i_export[1] += 1
        if mod(i_export[1], params[:export_frequency]) == 0
            @info "Exporting file: $(nc_q.filename)"
            append_data(nc_q, q, t[1])
            append_data(nc_aux, aux, t[1])
        end
    end
    extrap_0th_order!(q, (:θ_liq, :q_tot), grid, gm)
    extrap_0th_order!(aux, :T, grid, gm)

    append_data(nc_q, q, t[1])
    append_data(nc_aux, aux, t[1])

    return (grid, q, aux, params)
end
