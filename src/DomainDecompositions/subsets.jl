#### Domain subsets

export DomainSubSet
export get_param

"""
    DomainSubSet{GM,EN,UD}

Subset of domains
 - `gm` bool indicating to include grid-mean domain
 - `en` bool indicating to include environment sub-domain
 - `ud` bool indicating to include updraft sub-domains
"""
struct DomainSubSet{GM,EN,UD}
  DomainSubSet(;gm::GM=false,en::EN=false,ud::UD=false) where {GM<:Bool,EN<:Bool,UD<:Bool} =
    new{GridMean{gm},Environment{en},Updraft{ud}}()
end

gridmean(   ::DomainSubSet{GM,EN,UD}) where {GM,EN,UD} = n_subdomains(GM)
environment(::DomainSubSet{GM,EN,UD}) where {GM,EN,UD} = n_subdomains(EN)
updraft(    ::DomainSubSet{GM,EN,UD}) where {GM,EN,UD} = n_subdomains(UD)

gridmean(   ::DomainSet{GM,EN,UD}, domain_subset::DomainSubSet) where {GM,EN,UD} =
  gridmean(domain_subset)    ? n_subdomains(GM) : 0
environment(::DomainSet{GM,EN,UD}, domain_subset::DomainSubSet) where {GM,EN,UD} =
  environment(domain_subset) ? n_subdomains(EN) : 0
updraft(    ::DomainSet{GM,EN,UD}, domain_subset::DomainSubSet) where {GM,EN,UD} =
  updraft(domain_subset)     ? n_subdomains(UD) : 0

n_subdomains(domain_set::DomainSet, domain_subset::DomainSubSet) =
  (gridmean(domain_set, domain_subset), environment(domain_set, domain_subset), updraft(domain_set, domain_subset))
Base.sum( domain_set::DomainSet, domain_subset::DomainSubSet) =
  gridmean(domain_set, domain_subset) + environment(domain_set, domain_subset) + updraft(domain_set, domain_subset)
