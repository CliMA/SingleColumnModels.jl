#### ProcessResults

# TODO: Migrate to using NetCDF IO
using NCDatasets

function export_data(q, aux, grid, output_dir, params)
    if params[:export_data]
        nc_q = NetCDFWriter(joinpath(output_dir, "q"))
        nc_aux = NetCDFWriter(joinpath(output_dir, "aux"))
        export_state(nc_q, grid, q)
        export_state(nc_aux, grid, q)
    end
end
