"""
    DomainDecompositions

Provides a domain decomposition into
 - `grid_mean`
 - `environment`
 - `updrafts`
"""
module DomainDecompositions

export DomainDecomposition
export gridmean, environment, updraft
export over_sub_domains, allcombinations
export var_names, var_string, var_suffix

include("types.jl")
include("subsets.jl")
include("indexes.jl")
include("variable_mapper.jl")

"""
    DomainDecomposition(vars, domain_set::DomainSet) where {T}

Return a domain decomposition, given a tuple of tuples of variable
names and the number of their subdomains.
"""
function DomainDecomposition(vars, domain_set::DomainSet) where {T}
    n_subdomains = sum(domain_set)
    n_vars = sum([sum(domain_set, dss) for (v, dss) in vars])
    var_mapper, dss_per_var, var_names = get_var_mapper(vars, domain_set)
    idx = DomainIdx(domain_set)
    sd_unmapped = get_sd_unmapped(vars, idx, domain_set)

    idx_ss_per_var =
        Dict(ϕ => DomainIdx(domain_set, dss_per_var[ϕ]) for ϕ in var_names)
    a_map = Dict(ϕ => get_sv_a_map(idx, idx_ss_per_var[ϕ]) for ϕ in var_names)

    return DomainDecomposition{typeof(var_mapper)}(
        idx,
        var_mapper,
        idx_ss_per_var,
        a_map,
        sd_unmapped,
        domain_set,
        dss_per_var,
    )
end

gridmean(domain_decomp::DomainDecomposition) = gridmean(domain_decomp.idx)
environment(domain_decomp::DomainDecomposition) = environment(domain_decomp.idx)
updraft(domain_decomp::DomainDecomposition) = updraft(domain_decomp.idx)
alldomains(domain_decomp::DomainDecomposition) = alldomains(domain_decomp.idx)

"""
    over_sub_domains(domain_decomp::DomainDecomposition, ϕ::Symbol)

Get list of indexes over all subdomains for variable `ϕ`.
"""
function over_sub_domains(domain_decomp::DomainDecomposition, ϕ::Symbol)
    return [x for x in domain_decomp.sd_unmapped[ϕ] if !(x == 0)]
end

"""
    var_names(domain_decomp::DomainDecomposition)

Get tuple of variable names.
"""
var_names(domain_decomp::DomainDecomposition{VM}) where {VM} = fieldnames(VM)

function Base.show(io::IO, domain_decomp::DomainDecomposition)
    println(io, "\n----------------------- Domain Decomposition")
    println(io, "DomainIdx")
    println(io, domain_decomp.idx)
    println(io, "var_mapper")
    println(io, domain_decomp.var_mapper)
    println(io, "\n-----------------------")
end

DomainIdx(domain_decomp::DomainDecomposition) = domain_decomp.idx

allcombinations(domain_decomp::DomainDecomposition) =
    allcombinations(DomainIdx(domain_decomp))

end
