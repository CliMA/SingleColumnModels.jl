#### Domain decomposition
export DomainSet

import Base
abstract type Domain{N} end
struct GridMean{N} <: Domain{N} end
struct Environment{N} <: Domain{N} end
struct Updraft{N} <: Domain{N} end

"""
    n_subdomains(::Type{Domain})

The number of domains in domain `Domain`
"""
n_subdomains(::Type{GridMean{N}}) where N = N
n_subdomains(::Type{Environment{N}}) where N = N
n_subdomains(::Type{Updraft{N}}) where N = N

"""
    DomainSet{GM,EN,UD}

Decomposition of number of domains
 - `gm` number of grid-mean domain
 - `en` number of environment sub-domain
 - `ud` number of updraft sub-domains
"""
struct DomainSet{GM,EN,UD}
  function DomainSet(;gm::GM=0,en::EN=0,ud::UD=0) where {GM<:Int,EN<:Int,UD<:Int}
    @assert 0<=gm<=1
    @assert 0<=en<=1
    @assert 0<=ud
    return new{GridMean{gm},Environment{en},Updraft{ud}}()
  end
end

Base.sum(   ::DomainSet{GM,EN,UD}) where {GM,EN,UD} = n_subdomains(GM)+n_subdomains(EN)+n_subdomains(UD)
n_subdomains(  ::DomainSet{GM,EN,UD}) where {GM,EN,UD} = (n_subdomains(GM),n_subdomains(EN),n_subdomains(UD))

"""
    DomainIdx{GM,EN,UD}

Decomposition of domain indexes, including indexes for
 - `i_gm` grid-mean domain
 - `i_en` environment sub-domain
 - `i_ud` updraft sub-domains

Note that the index ordering logic is defined here.
"""
struct DomainIdx{GM,EN,UD} end

"""
    DomainDecomposition{VM}

A domain decomposition into
 - `grid_mean`
 - `environment`
 - `updrafts`
"""
struct DomainDecomposition{VM}
  """
  A `DomainIdx`, containing the complete set of sub-domain indexes
  """
  idx::DomainIdx
  """
  A `NamedTuple` variable mapper, containing indexes per variable
  """
  var_mapper::VM
  """
  A `Dict` index per sub-domain per variable
  """
  idx_ss_per_var::Dict
  """
  A mapper, containing mapping from global sub-domain indexes to sub-domain indexes per variable
  """
  a_map::Dict
  """
  A `NamedTuple` sub-domain mapper, containing indexes per variable which are sub-domain specific
  """
  sd_unmapped::Dict
  """
  A `DomainSet`, which specifies the complete set of `DomainSet` across all variables
  """
  domain_set::DomainSet
  """
  A `DomainSet` per variable
  """
  dss_per_var::NamedTuple
end
