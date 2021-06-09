#### ProcessResults

# TODO: Migrate to using NetCDF IO
using NCDatasets

function export_data(q, tmp, grid, output_dir, params)
    if params[:export_data]
        nc_q = NetCDFWriter(joinpath(output_dir, "q"))
        nc_tmp = NetCDFWriter(joinpath(output_dir, "tmp"))
        export_state(nc_q, grid, q)
        export_state(nc_tmp, grid, q)
    end
end
