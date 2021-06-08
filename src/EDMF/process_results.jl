#### ProcessResults

# TODO: Migrate to using NetCDF IO
using NCDatasets

function export_data(q, tmp, grid, dir_tree, params)
    if params[:export_data]
        nc_q = NetCDFWriter(joinpath(dir_tree.output, "q"))
        nc_tmp = NetCDFWriter(joinpath(dir_tree.output, "tmp"))
        export_state(nc_q, grid, q)
        export_state(nc_tmp, grid, q)
    end
end
