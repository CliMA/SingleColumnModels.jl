#### TurbConvs

export TurbConv
export get_ϕ_ψ

struct TurbConv{G, SVQ, SVT, TD, OD}
    grid::G
    q::SVQ
    q_new::SVQ
    tendencies::SVQ
    aux::SVT
    tri_diag::TD
    aux_O2::Dict
    output_dir::OD
end

"""
    UniformGridParams{FT,IT}
Grid parameters used to initialize a grid.
"""
Base.@kwdef struct UniformGridParams{FT, IT}
    z_min::FT
    z_max::FT
    n_elems::IT
end

import ..FiniteDifferenceGrids: UniformGrid

UniformGrid(gp::UniformGridParams) = UniformGrid(gp.z_min, gp.z_max, gp.n_elems)


"""
    get_ϕ_ψ(ϕ)

Convenience function to get individual
variables names from co-variance.
"""
function get_ϕ_ψ(ϕ)
    if ϕ == :cv_q_tot
        _ϕ = :q_tot
        _ψ = :q_tot
    elseif ϕ == :cv_θ_liq
        _ϕ = :θ_liq
        _ψ = :θ_liq
    elseif ϕ == :cv_θ_liq_q_tot
        _ϕ = :q_tot
        _ψ = :θ_liq
    end
    return _ϕ, _ψ
end

function TurbConv(params, case::Case, output_dir::String)
    @unpack N_subdomains = params
    n_ud = N_subdomains - 2

    grid = UniformGrid(params[:UniformGridParams])

    #! format: off

    domain_set = DomainSet(gm = 1, en = 1, ud = n_ud)

    unkowns = (
        (:a             , DomainSubSet(gm = true, en = true, ud = true)),
        (:w             , DomainSubSet(gm = true, en = true, ud = true)),
        (:q_tot         , DomainSubSet(gm = true, en = true, ud = true)),
        (:θ_liq         , DomainSubSet(gm = true, en = true, ud = true)),
        (:tke           , DomainSubSet(en = true)),
        (:u             , DomainSubSet(gm = true, en = true, ud = true)),
        (:v             , DomainSubSet(gm = true, en = true, ud = true)),
    )

    aux_vars = (
        (:ρ_0           , DomainSubSet(gm = true)),
        (:p_0           , DomainSubSet(gm = true)),
        (:T             , DomainSubSet(gm = true, en = true, ud = true)),
        (:q_liq         , DomainSubSet(gm = true, en = true, ud = true)),
        (:HVSD_a        , DomainSubSet(ud = true)),
        (:HVSD_w        , DomainSubSet(ud = true)),
        (:buoy          , DomainSubSet(gm = true, en = true, ud = true)),
        (:nh_press      , DomainSubSet(ud = true)),
        (:∇buoyancy     , DomainSubSet(gm = true)),
        (:δ_model       , DomainSubSet(ud = true)),
        (:ε_model       , DomainSubSet(ud = true)),
        (:l_mix         , DomainSubSet(gm = true)),
        (:K_m           , DomainSubSet(gm = true)),
        (:K_h           , DomainSubSet(gm = true)),
        (:α_0           , DomainSubSet(gm = true)),
        (:q_tot_dry     , DomainSubSet(gm = true)),
        (:θ_dry         , DomainSubSet(gm = true)),
        (:subsidence    , DomainSubSet(gm = true)),
        (:ug            , DomainSubSet(gm = true)),
        (:vg            , DomainSubSet(gm = true)),
        (:dTdt          , DomainSubSet(gm = true)),
        (:dqtdt         , DomainSubSet(gm = true)),
        (:t_cloudy      , DomainSubSet(gm = true)),
        (:q_vap_cloudy  , DomainSubSet(gm = true)),
        (:q_tot_cloudy  , DomainSubSet(gm = true)),
        (:θ_cloudy      , DomainSubSet(gm = true)),
        (:CF            , DomainSubSet(gm = true)),
        (:mf_θ_liq      , DomainSubSet(gm = true)),
        (:mf_q_tot      , DomainSubSet(gm = true)),
        (:mf_tend_θ_liq , DomainSubSet(gm = true)),
        (:mf_tend_q_tot , DomainSubSet(gm = true)),
        (:mf_aux        , DomainSubSet(ud = true)),
        (:θ_ρ           , DomainSubSet(gm = true)),
    )

    tri_diag_vars = (
        (:a             , DomainSubSet(gm = true)),
        (:b             , DomainSubSet(gm = true)),
        (:c             , DomainSubSet(gm = true)),
        (:f             , DomainSubSet(gm = true)),
        (:ρaK           , DomainSubSet(gm = true)),
    )

    q_2MO_vars = (
        (:entr_gain     , DomainSubSet(gm = true)),
        (:buoy          , DomainSubSet(gm = true)),
        (:press         , DomainSubSet(gm = true)),
        (:shear         , DomainSubSet(gm = true)),
        (:interdomain   , DomainSubSet(gm = true)),
    )
    #! format: on

    q = StateVec(unkowns, grid, domain_set)
    aux = StateVec(aux_vars, grid, domain_set)
    tri_diag = StateVec(tri_diag_vars, grid, domain_set)
    q_new = deepcopy(q)
    tendencies = deepcopy(q)
    aux_O2 = Dict()
    aux_O2[:tke] = StateVec(q_2MO_vars, grid, domain_set)

    turb_conv =
        TurbConv(grid, q, q_new, tendencies, aux, tri_diag, aux_O2, output_dir)
end
