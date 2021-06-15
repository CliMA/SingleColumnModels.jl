"""
    StateVecs

Provides a state vector that is includes information
of the domain decomposition into a grid mean,
environment and updrafts.
"""
module StateVecs

using ClimateMachineCore.Spaces

# Forwarded methods
import ..DomainDecompositions:
    var_names,
    DomainIdx,
    allcombinations,
    gridmean,
    environment,
    updraft,
    over_sub_domains,
    alldomains,
    get_i_var,
    var_string,
    var_suffix
using ..DomainDecompositions

export StateVec, over_sub_domains, Cut, Dual, allcombinations
export var_names, var_string, var_suffix
export assign!, assign_real!, assign_ghost!, extrap!, extrap_0th_order!
export compare

struct FieldsPerElement
    vals::Vector
end

"""
    StateVec{VM}

A state vector containing the number of subdomains,
`n_subdomains`, a `NamedTuple` variable mapper, a tuple
of the variable names, and a vector of vectors, containing
the values for all of the variables.
"""
struct StateVec{VM}
    """
    A `DomainDecomposition`, containing the map between subdomains variables and the state vector index
    """
    domain_decomp::DomainDecomposition{VM}
    """
    A vector of vectors containing all values for all variables.
    """
    fields::Vector{FieldsPerElement}
end

"""
    StateVec(vars::Tuple{Vararg{Tuple{Symbol, DD}}}, grid::Grid{T}) where {DD, T}

Return a state vector, given a tuple of tuples of variable
names and the number of their subdomains.
"""
function StateVec(
    vars,
    ::Type{FT},
    n_points::Int,
    domain_set::DomainSet,
) where {FT}
    domain_decomp = DomainDecomposition(vars, domain_set)
    n_vars = sum([sum(domain_set, dss) for (v, dss) in vars])
    all_vars = Any[FT(0) for v in 1:n_vars]
    fields = [FieldsPerElement(deepcopy(all_vars)) for k in 1:n_points]
    return StateVec{typeof(domain_decomp.var_mapper)}(domain_decomp, fields)
end


function Base.show(io::IO, sv::StateVec)
    println(io, "\n----------------------- State Vector")
    println(io, sv.domain_decomp)
    vn = var_names(sv)
    headers = [
        length(over_sub_domains(sv, ϕ)) == 1 ? string(ϕ) :
            string(ϕ) * "_" * string(i) for ϕ in vn
        for i in over_sub_domains(sv, ϕ)
    ]
    n_vars = length(headers)
    n_elem = length(sv.fields)
    data = reshape(
        [
            sv[ϕ, k, i] for ϕ in vn for i in over_sub_domains(sv, ϕ)
            for k in 1:n_elem
        ],
        n_elem,
        n_vars,
    )
    Base.print_matrix(io, [reshape(Text.(headers), 1, n_vars); data])
    println(io, "\n-----------------------")
end

# Overload
var_names(sv::StateVec) = var_names(sv.domain_decomp)
alldomains(sv::StateVec) = alldomains(sv.domain_decomp)
DomainIdx(sv::StateVec) = DomainIdx(sv.domain_decomp)
allcombinations(sv::StateVec) = allcombinations(sv.domain_decomp)
gridmean(sv::StateVec) = gridmean(sv.domain_decomp)
environment(sv::StateVec) = environment(sv.domain_decomp)
updraft(sv::StateVec) = updraft(sv.domain_decomp)
get_i_var(sv::StateVec, ϕ::Symbol, i) = get_i_var(sv.domain_decomp, ϕ, i)

"""
    over_sub_domains(sv::StateVec, ϕ::Symbol)

Get list of indexes over all subdomains for variable `ϕ`.
"""
over_sub_domains(sv::StateVec, ϕ::Symbol) =
    over_sub_domains(sv.domain_decomp, ϕ)

function Base.getindex(sv::StateVec, ϕ::Symbol, k, i = 0)
    i == 0 && (i = gridmean(sv))
    i_sv = DomainDecompositions.get_i_state_vec(sv.domain_decomp, ϕ, i)
    @inbounds sv.fields[k].vals[i_sv]
end
function Base.setindex!(sv::StateVec, val, ϕ::Symbol, k, i = 0)
    i == 0 && (i = gridmean(sv))
    @inbounds i_sv =
        DomainDecompositions.get_i_state_vec(sv.domain_decomp, ϕ, i)
    @inbounds sv.fields[k].vals[i_sv] = val
end

abstract type AbstractCut{I} end

"""
    Cut{I} <: AbstractCut{I}

A Cut struct used to slice the state
vector along the grid-element dimension.
This is used as an API to pass Cuts
into local derivative/interpolation routines.
"""
struct Cut{I} <: AbstractCut{I}
    e::I
end
function Base.getindex(sv::StateVec, ϕ::Symbol, cut::Cut, i = 0)
    @inbounds [sv[ϕ, k, i] for k in (cut.e - 1):(cut.e + 1)]
end

"""
    Dual{I} <: AbstractCut{I}

A Dual struct used to slice the state
vector along the grid-element dimension.
In addition, this interpolates fields
[(f[k-1]+f[k])/2, (f[k]+f[k+1])/2]
"""
struct Dual{I} <: AbstractCut{I}
    e::I
end
function Base.getindex(sv::StateVec, ϕ::Symbol, dual::Dual, i = 0)
    @inbounds [
        (sv[ϕ, k, i] + sv[ϕ, k + 1, i]) / 2 for k in (dual.e - 1):(dual.e)
    ]
end

Base.isnan(sv::StateVec) = any([any(isnan.(fpe.vals)) for fpe in sv.fields])
Base.isinf(sv::StateVec) = any([any(isinf.(fpe.vals)) for fpe in sv.fields])


"""
    var_suffix(sv::StateVec, ϕ::Symbol, i=0)

A suffix string for variable `ϕ` indicating its sub-domain.
"""
var_suffix(sv::StateVec, ϕ::Symbol, i = 0) = var_suffix(sv.domain_decomp, ϕ, i)

"""
    var_string(sv::StateVec, ϕ::Symbol, i=0)

A string of the variable `ϕ` with the appropriate suffix indicating its sub-domain.
"""
var_string(sv::StateVec, ϕ::Symbol, i = 0) = var_string(sv.domain_decomp, ϕ, i)

"""
    assign!(sv::StateVec, var_names, grid::Grid, val, i=1)

Assign value `val` to variable `ϕ` for all ghost points.
"""
function assign!(sv::StateVec, var_names, grid::Grid{FT}, val::FT) where {FT}
    gm, en, ud, sd, al = allcombinations(sv)
    !(var_names isa Tuple) && (var_names = (var_names,))
    @inbounds for k in over_elems(grid),
        ϕ in var_names,
        i in over_sub_domains(sv, ϕ)

        sv[ϕ, k, i] = val
    end
end

"""
    assign!(sv::StateVec, grid::Grid, val)

Assign value `val` to all variables in state vector.
"""
function assign!(sv::StateVec, grid::Grid{FT}, val::FT) where {FT}
    gm, en, ud, sd, al = allcombinations(sv)
    @inbounds for k in over_elems(grid),
        ϕ in var_names(sv),
        i in over_sub_domains(sv, ϕ)

        sv[ϕ, k, i] = val
    end
end

"""
    assign_real!(sv::StateVec, ϕ::Symbol, grid::Grid, vec_real, i=0)

Assign the real elements to the given vector `vec_real`.
"""
function assign_real!(sv::StateVec, ϕ::Symbol, grid::Grid, vec_real, i = 0)
    i == 0 && (i = gridmean(sv))
    k_real = 1
    @assert length(vec_real) == length(over_elems_real(grid))
    @inbounds for k in over_elems_real(grid)
        sv[ϕ, k, i] = vec_real[k_real]
        k_real += 1
    end
end

"""
    assign_ghost!(sv::StateVec, ϕ::Symbol, grid::Grid, val, i=0)

Assign value `val` to variable `ϕ` for all ghost points.
"""
function assign_ghost!(sv::StateVec, ϕ::Symbol, grid::Grid, val, i = 0)
    i == 0 && (i = gridmean(sv))
    @inbounds for k in over_elems_ghost(grid)
        sv[ϕ, k, i] = val
    end
end

"""
    assign_ghost!(sv::StateVec, ϕ::Symbol, grid::Grid, val, b::ZBoundary, i=1)

Assign value `val` to variable `ϕ` for all ghost points.
"""
function assign_ghost!(
    sv::StateVec,
    ϕ::Symbol,
    grid::Grid,
    val,
    b::ZBoundary,
    i = 0,
)
    i == 0 && (i = gridmean(sv))
    @inbounds for k in over_elems_ghost(grid, b)
        sv[ϕ, k, i] = val
    end
end

"""
    extrap!(sv::StateVec, ϕ::Symbol, grid::Grid, i=1)

Extrapolate variable `ϕ` to the first ghost point.
"""
function extrap!(sv::StateVec, ϕ::Symbol, grid::Grid, i = 0)
    i == 0 && (i = gridmean(sv))
    @inbounds for b in (Zmin(), Zmax())
        kg, ki, kii = boundary_points(grid, b)
        sv[ϕ, kg, i] = 2 * sv[ϕ, ki, i] - sv[ϕ, kii, i]
    end
end

"""
    extrap!(sv::StateVec, ϕ::Symbol, grid::Grid, ::ZBoundary, i=1)

Extrapolate variable `ϕ` to the first ghost point.
"""
function extrap!(sv::StateVec, ϕ::Symbol, grid::Grid, b::ZBoundary, i = 0)
    i == 0 && (i = gridmean(sv))
    kg, ki, kii = boundary_points(grid, b)
    @inbounds sv[ϕ, kg, i] = 2 * sv[ϕ, ki, i] - sv[ϕ, kii, i]
end

"""
    extrap!(sv::StateVec, ϕ::Symbol, grid::Grid, dual_val, b::ZBoundary, i=1)

Extrapolate variable `ϕ` to the first ghost point.
"""
function extrap!(
    sv::StateVec,
    ϕ::Symbol,
    grid::Grid,
    dual_val,
    b::ZBoundary,
    i = 0,
)
    i == 0 && (i = gridmean(sv))
    ki = first_interior(grid, b)
    kg = first_ghost(grid, b)
    @inbounds sv[ϕ, kg, i] = 2 * dual_val - sv[ϕ, ki, i]
end

"""
    extrap_0th_order!(sv::StateVec, ϕ::Symbol, grid::Grid, i=1)

Extrapolate variable `ϕ` to the first ghost point using zeroth order approximation.
"""
function extrap_0th_order!(sv::StateVec, var_names, grid::Grid, i = 0)
    !(var_names isa Tuple) && (var_names = (var_names,))
    i == 0 && (i = gridmean(sv))
    @inbounds for ϕ in var_names
        @inbounds for b in (Zmin(), Zmax())
            kg, ki, kii = boundary_points(grid, b)
            sv[ϕ, kg, i] = sv[ϕ, ki, i]
        end
    end
end

include("boundary_conditions.jl")
include("export_funcs.jl")
include("domain_decomposition_funcs.jl")

end
