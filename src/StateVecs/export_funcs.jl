#### ExportFuncs

export NetCDFWriter
export export_state
export append_data
export init_data
export full_name

using NCDatasets
using OrderedCollections

#####
##### Diagnostic variable map
#####

struct DiagnosticVariable
    name::String
    attrib::OrderedDict
    DiagnosticVariable(name::String, attrib::OrderedDict = OrderedDict()) =
        new(name, attrib)
end

function var_attrib(
    units::String = "[units]",
    long_name::String = "[long_name_default]",
    standard_name::String = "[standard_name_default]",
    fill_value::Union{Nothing, Any} = nothing,
)
    attrib = OrderedDict(
        "units" => units,
        "long_name" => long_name,
        "standard_name" => standard_name,
    )
    if !isnothing(fill_value)
        attrib["_FillValue"] = fill_value
    end
    return attrib
end

#! format: off
const GLOBAL_VARIABLE_MAP = OrderedDict{String, DiagnosticVariable}(
    "w"      => DiagnosticVariable("w"     , var_attrib("m s^-1", "vertical wind", "upward_air_velocity")),
    "ρ"      => DiagnosticVariable("ρ"     , var_attrib("kg m^-3", "air density", "air_density")),
    "a_gm"   => DiagnosticVariable("ρ"     , var_attrib("kg m^-3", "air density", "air_density")),
    "ρ_0_gm" => DiagnosticVariable("ρ_0_gm", var_attrib("units", "ρ_0_gm", "ρ_0_gm")),
    "a_gm"   => DiagnosticVariable("a_gm"  , var_attrib("units", "a_gm", "a_gm")),
    "a_en"   => DiagnosticVariable("a_en"  , var_attrib("units", "a_en", "a_en")),
    "a_ud_1" => DiagnosticVariable("a_ud_1", var_attrib("units", "a_ud_1", "a_ud_1")),
    "w_gm"   => DiagnosticVariable("w_gm"  , var_attrib("units", "w_gm", "w_gm")),
    "w_en"   => DiagnosticVariable("w_en"  , var_attrib("units", "w_en", "w_en")),
    "w_ud_1" => DiagnosticVariable("w_ud_1", var_attrib("units", "w_ud_1", "w_ud_1")),
    "ϕ_gm"   => DiagnosticVariable("ϕ_gm"  , var_attrib("units", "ϕ_gm", "ϕ_gm")),
    "ϕ_en"   => DiagnosticVariable("ϕ_en"  , var_attrib("units", "ϕ_en", "ϕ_en")),
    "ϕ_ud_1" => DiagnosticVariable("ϕ_ud_1", var_attrib("units", "ϕ_ud_1", "ϕ_ud_1")),
    "ψ_gm"   => DiagnosticVariable("ψ_gm"  , var_attrib("units", "ψ_gm", "ψ_gm")),
    "ψ_en"   => DiagnosticVariable("ψ_en"  , var_attrib("units", "ψ_en", "ψ_en")),
    "ψ_ud_1" => DiagnosticVariable("ψ_ud_1", var_attrib("units", "ψ_ud_1", "ψ_ud_1")),
    "ρ_0_gm" => DiagnosticVariable("ρ_0_gm", var_attrib("units", "ρ_0_gm", "ρ_0_gm"))
)
#! format: on

#####
##### NetCDF writer
#####

abstract type AbstractWriter end

struct NetCDFWriter{M} <: AbstractWriter
    filename::String
    meta::M
end
NetCDFWriter(filename::String) = NetCDFWriter{Nothing}(filename, nothing)

full_name(writer::NetCDFWriter) = writer.filename * ".nc"

function get_attrib(var_name)
    if haskey(GLOBAL_VARIABLE_MAP, var_name)
        var = GLOBAL_VARIABLE_MAP[var_name]
        return var.attrib
    else
        return var_attrib()
    end
end

#####
##### Implementation
#####

function export_state(
    nc::NetCDFWriter,
    grid,
    sv::StateVec,
    simtime = 0;
    kwargs...,
)
    nc = init_data(nc, grid, sv; kwargs...)
    append_data(nc, sv, simtime)
end

function init_data(nc::NetCDFWriter, grid, sv::StateVec; kwargs...)
    # TODO: allow exporting onto cell faces
    vars = OrderedDict()
    include_ghost = true
    domain_range = include_ghost ? over_elems(grid) : over_elems_real(grid)
    z = grid.zc[domain_range]
    dims = OrderedDict("z" => (z, OrderedDict()))

    # For collecting var_names and adding to GLOBAL_VARIABLE_MAP:
    # @inbounds for name_id in var_names(sv)
    #     @inbounds for i in over_sub_domains(sv, name_id)
    #         var_name = var_string(sv, name_id, i)
    #         @show var_name
    #     end
    # end

    FT = eltype(grid)

    @inbounds for name_id in var_names(sv)
        @inbounds for i in over_sub_domains(sv, name_id)
            var_name = var_string(sv, name_id, i)
            attrib = get_attrib(var_name)
            vars[var_name] = (tuple(collect(keys(dims))...), FT, attrib)
        end
    end

    init_data(nc, dims, vars; kwargs...)
    return NetCDFWriter(nc.filename, domain_range)
end

function init_data(nc::NetCDFWriter, dims, vars; overwrite = true)
    nm = full_name(nc)
    if isfile(nm)
        if overwrite
            @warn "$(nm) exists and will be overwritten."
        else
            error("$(nm) exists and overwriting is forbidden.")
        end
    end
    Dataset(nm, "c") do ds
        # define spatial and time dimensions
        for (dn, (dv, da)) in dims
            defDim(ds, dn, length(dv))
        end
        defDim(ds, "time", Inf) # Inf sets UNLIMITED dimension

        # include dimensions as variables
        for (dn, (dv, da)) in dims
            defVar(ds, dn, dv, (dn,), attrib = da)
        end
        defVar(
            ds,
            "time",
            Float64,
            ("time",),
            attrib = OrderedDict(
                "units" => "seconds since 1900-01-01 00:00:00",
                "long_name" => "time",
            ),
        )

        # define variables
        for (vn, (vd, vt, va)) in vars
            defVar(ds, vn, vt, (vd..., "time"), attrib = va)
        end
    end
    return nothing
end

function append_data(nc::NetCDFWriter, sv::StateVec, simtime = 0)
    varvals = OrderedDict()
    domain_range = nc.meta
    @inbounds for name_id in var_names(sv)
        @inbounds for i in over_sub_domains(sv, name_id)
            var_name = var_string(sv, name_id, i)
            varvals[var_name] = [sv[name_id, k, i] for k in domain_range]
        end
    end
    append_data(nc, varvals, simtime)
end

function append_data(nc::NetCDFWriter, varvals, simtime)
    Dataset(full_name(nc), "a") do ds
        timevar = ds["time"]
        t = length(timevar) + 1
        timevar[t] = simtime
        for (vn, vv) in varvals
            dsvar = ds[vn]
            dsvar[ntuple(_ -> Colon(), ndims(vv))..., t] = vv
        end
    end
    return nothing
end
