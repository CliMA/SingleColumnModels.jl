#### RefState

using OrdinaryDiffEq

"""
    init_ref_state!(aux::StateVec,
                    grid::Grid,
                    params,
                    output_dir::String)

Initializes the reference state variables:
  - `p_0` pressure
  - `ρ_0` density
  - `α_0` specific volume
assuming a hydrostatic balance.
FIXME: add reference
"""
function init_ref_state!(
    aux::StateVec,
    grid::Grid{FT},
    params,
    output_dir::String,
) where {FT}
    @unpack param_set, SurfaceModel = params
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
    assign_real!(aux, :p_0, grid, p_0)
    apply_Neumann!(aux, :p_0, grid, 0.0, Zmin())
    apply_Neumann!(aux, :p_0, grid, 0.0, Zmax())

    @inbounds for k in over_elems(grid)
        ts = PhaseEquil_pTq(param_set, aux[:p_0, k], T_g, q_tot_g)
        q_pt = PhasePartition(ts)
        T = air_temperature(ts)
        aux[:ρ_0, k] = air_density(param_set, T, aux[:p_0, k], q_pt)
        aux[:α_0, k] = 1 / aux[:ρ_0, k]
    end
    extrap!(aux, :ρ_0, grid)
    extrap!(aux, :α_0, grid)
    extrap!(aux, :p_0, grid)

    nc = NetCDFWriter(joinpath(output_dir, "aux_ref_state"))
    export_state(nc, grid, aux)
end
