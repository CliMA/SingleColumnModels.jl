#### RefState

using OrdinaryDiffEq

"""
    init_ref_state!(tmp::StateVec,
                    grid::Grid,
                    params,
                    dir_tree::DirTree)

Initializes the reference state variables:
  - `p_0` pressure
  - `ρ_0` density
  - `α_0` specific volume
assuming a hydrostatic balance.
FIXME: add reference
"""
function init_ref_state!(
    tmp::StateVec,
    grid::Grid{FT},
    params,
    dir_tree::DirTree,
) where {FT}
    @unpack params param_set SurfaceModel
    T_g = SurfaceModel.T
    q_tot_g = SurfaceModel.q_tot
    P_g = SurfaceModel.P

    q_pt_g = PhasePartition(q_tot_g)
    θ_liq_ice_g =
        TD.liquid_ice_pottemp_given_pressure(param_set, T_g, P_g, q_pt_g)
    logp = log(P_g)

    function tendencies(p, u, z)
        expp = exp(p)
        ρ = air_density(param_set, T_g, expp, q_pt_g)
        ts = LiquidIcePotTempSHumEquil_old(
            param_set,
            θ_liq_ice_g,
            q_tot_g,
            ρ,
            expp,
        )
        R_m = gas_constant_air(ts)
        T = air_temperature(ts)
        return -FT(grav(param_set)) / (T * R_m)
    end

    z_span = (grid.zn_min, grid.zn_max)
    prob = ODEProblem(tendencies, logp, z_span)
    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12)
    p_0 = [exp(sol(grid.zc[k])) for k in over_elems_real(grid)]
    assign_real!(tmp, :p_0, grid, p_0)
    apply_Neumann!(tmp, :p_0, grid, 0.0, Zmin())
    apply_Neumann!(tmp, :p_0, grid, 0.0, Zmax())

    @inbounds for k in over_elems(grid)
        ts = PhaseEquil_pTq(param_set, tmp[:p_0, k], T_g, q_tot_g)
        q_pt = PhasePartition(ts)
        T = air_temperature(ts)
        tmp[:ρ_0, k] = air_density(param_set, T, tmp[:p_0, k], q_pt)
        tmp[:α_0, k] = 1 / tmp[:ρ_0, k]
    end
    extrap!(tmp, :ρ_0, grid)
    extrap!(tmp, :α_0, grid)
    extrap!(tmp, :p_0, grid)

    plot_state(tmp, grid, dir_tree[:initial_conditions], :p_0)
    plot_state(tmp, grid, dir_tree[:initial_conditions], :ρ_0)
    plot_state(tmp, grid, dir_tree[:initial_conditions], :α_0)

end
