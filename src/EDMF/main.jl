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

    gm, en, ud, sd, al = allcombinations(q)

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

    nc_q = NetCDFWriter(joinpath(output_dir, "prog_vs_time"))
    nc_q = init_data(nc_q, grid, q)
    nc_aux = NetCDFWriter(joinpath(output_dir, "aux_vs_time"))
    nc_aux = init_data(nc_aux, grid, aux)

    @show joinpath(pwd(), full_name(nc_q))

    while t[1] < t_end
        assign!(q_tendencies, (:u, :v, :q_tot, :θ_liq), grid, 0.0)

        update_aux!(grid, q, aux, aux_O2, params, case)

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

        i_export[1] += 1
        if mod(i_export[1], params[:export_frequency]) == 0
            @info "Exporting file: $(nc_q.filename)"
            append_data(nc_q, q, t[1])
            append_data(nc_aux, aux, t[1])
        end
    end

    append_data(nc_q, q, t[1])
    append_data(nc_aux, aux, t[1])

    return (grid, q, aux, params)
end
