using Plots
using Statistics

function plot_profiles(ds_filename::String, output_dir; kwargs...)
    Dataset(ds_filename) do ds
        plot_profiles(ds, output_dir; kwargs...)
    end
end

function plot_profiles(
    ds,
    output_dir;
    skip_fields = [],
    fields_only = nothing,
    profiles_only = false,
)
    mkpath(output_dir)
    z = ds["z"][:]
    t = ds["time"][:]
    for k in keys(ds)
        k == "z" && continue
        k == "time" && continue
        k in skip_fields && continue
        if fields_only ≠ nothing
            k in fields_only || continue
        end

        vals = ds[k][:, :]
        file_ext = ".png"
        figname = "profile_$k$file_ext"

        plot(mean(vals, dims = 2), z; xlabel = k, ylabel = "height (m)")
        savefig(joinpath(output_dir, figname))

        profiles_only || continue
        figname = "contour_$k$file_ext"
        contourf(t, z, vals; xlabel = k, ylabel = "height (m)", c = :viridis)
        savefig(joinpath(output_dir, figname))
    end
end
