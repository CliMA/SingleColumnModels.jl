#### EDMF Functions

bound(x, x_bounds) = min(max(x, x_bounds[1]), x_bounds[2])

inside_bounds(x, x_bounds) = x > x_bounds[1] && x < x_bounds[2]


stable(obukhov_length::FT) where {FT} = obukhov_length > 0
unstable(obukhov_length::FT) where {FT} = obukhov_length < 0

"""
    StabilityDependentParam

A parameter that has two values:
 - `stable` for stable conditions
 - `unstable` for unstable conditions
"""
struct StabilityDependentParam{FT}
    stable::FT
    unstable::FT
end

function (sdp::StabilityDependentParam{FT})(obukhov_length::FT) where {FT}
    return unstable(obukhov_length) ? sdp.unstable : sdp.stable
end

function compute_tendencies_en_O2!(
    grid::Grid{FT},
    q_tendencies,
    aux_O2,
    cv,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q_tendencies)
    k_1 = first_interior(grid, Zmin())
    @inbounds for k in over_elems_real(grid)
        q_tendencies[cv, k, en] =
            aux_O2[cv][:press, k] +
            aux_O2[cv][:buoy, k] +
            aux_O2[cv][:shear, k] +
            aux_O2[cv][:entr_gain, k]
        # aux_O2[cv][:rain_src, k]
    end
    q_tendencies[cv, k_1, en] = FT(0)
end

function compute_tendencies_gm_scalars!(grid, q_tendencies, q, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack SurfaceModel = params
    k_1 = first_interior(grid, Zmin())
    Δzi = grid.Δzi
    α_1 = aux[:α_0, k_1]
    ae_1 = q[:a, k_1, en]
    @inbounds for k in over_elems(grid)
        q_tendencies[:q_tot, k, gm] += aux[:mf_tend_q_tot, k]
        q_tendencies[:θ_liq, k, gm] += aux[:mf_tend_θ_liq, k]
    end

    q_tendencies[:q_tot, k_1, gm] += SurfaceModel.ρq_tot_flux * Δzi * α_1 / ae_1
    q_tendencies[:θ_liq, k_1, gm] += SurfaceModel.ρθ_liq_flux * Δzi * α_1 / ae_1
end

function compute_tendencies_ud!(grid, q_tendencies, q, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            a_k = q[:a, k, i]
            α_0_kp = aux[:α_0, k]
            w_env = q[:w, k, en]
            ρ_k = aux[:ρ_0, k]
            w_i = q[:w, k, i]
            ε_model = aux[:ε_model, k, i]
            δ_model = aux[:δ_model, k, i]
            θ_liq_env = q[:θ_liq, k, en]
            q_tot_env = q[:q_tot, k, en]
            B_k = aux[:buoy, k, i]
            ρa_k = ρ_k * a_k
            ρaw_k = ρa_k * w_i

            a_cut = q[:a, Cut(k), i]
            θ_liq_cut = q[:θ_liq, Cut(k), i]
            q_tot_cut = q[:q_tot, Cut(k), i]
            w_cut = q[:w, Cut(k), i]
            ρ_cut = aux[:ρ_0, Cut(k)]

            ρaw_cut = ρ_cut .* a_cut .* w_cut
            ρawθ_liq_cut = ρaw_cut .* θ_liq_cut
            ρawq_tot_cut = ρaw_cut .* q_tot_cut
            ρaww_cut = ρ_cut .* a_cut .* w_cut .* w_cut

            tendencies = 0
            adv = -α_0_kp * advect(ρaw_cut, w_cut, grid)
            tendencies += adv
            ε_term = a_k * w_i * ε_model
            tendencies += ε_term
            δ_term = -a_k * w_i * δ_model
            tendencies += δ_term
            q_tendencies[:a, k, i] = tendencies

            adv = -advect(ρaww_cut, w_cut, grid)
            exch = ρaw_k * (-δ_model * w_i + ε_model * w_env)
            buoy = ρa_k * B_k
            nh_press = aux[:nh_press, k, i]

            tendencies = (adv + exch + buoy + nh_press)
            q_tendencies[:w, k, i] = tendencies

            tendencies_θ_liq = 0
            tendencies_q_tot = 0

            tendencies_θ_liq += -advect(ρawθ_liq_cut, w_cut, grid)
            tendencies_q_tot += -advect(ρawq_tot_cut, w_cut, grid)

            tendencies_θ_liq +=
                ρaw_k * (ε_model * θ_liq_env - δ_model * θ_liq_cut[2])
            tendencies_q_tot +=
                ρaw_k * (ε_model * q_tot_env - δ_model * q_tot_cut[2])

            q_tendencies[:θ_liq, k, i] = tendencies_θ_liq
            q_tendencies[:q_tot, k, i] = tendencies_q_tot
        end
    end
end

function saturation_adjustment_sd!(grid, q, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    @unpack param_set = params
    @inbounds for i in sd
        @inbounds for k in over_elems_real(grid)
            ts = ActiveThermoState(param_set, q, aux, k, i)
            q_liq = PhasePartition(ts).liq
            T = air_temperature(ts)
            aux[:T, k, i] = T
            aux[:q_liq, k, i] = q_liq
        end
    end
end

include(joinpath("Models", "entr_detr.jl"))
include(joinpath("Models", "buoyancy.jl"))
include(joinpath("Models", "mixing_length.jl"))
include(joinpath("Models", "microphysics.jl"))
include(joinpath("Models", "surface.jl"))
include(joinpath("Models", "pressure.jl"))
include(joinpath("Models", "eddy_diffusivity.jl"))

function top_of_updraft(grid::Grid, q::StateVec, params)
    @unpack w_bounds, a_bounds = params
    gm, en, ud, sd, al = allcombinations(q)

    z_star_a = zeros(length(ud))
    z_star_w = zeros(length(ud))
    k_2 = first_interior(grid, Zmax())
    @inbounds for i in ud
        z_star_a[i] = min([
            !inside_bounds(q[:a, k, i], a_bounds) ? grid.zc[k] :
                grid.zc[k_2 + 1] for k in over_elems_real(grid)
        ]...)
        z_star_w[i] = min([
            !inside_bounds(q[:w, k, i], w_bounds) ? grid.zc[k] :
                grid.zc[k_2 + 1] for k in over_elems_real(grid)
        ]...)
    end
    return z_star_a, z_star_w
end

function filter_scalars!(grid, q, aux, params)
    gm, en, ud, sd, al = allcombinations(q)
    z_star_a, z_star_w = top_of_updraft(grid, q, params)
    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            aux[:HVSD_a, k, i] = 1 - heaviside(grid.zc[k] - z_star_a[i], 1)
            aux[:HVSD_w, k, i] = 1 - heaviside(grid.zc[k] - z_star_w[i], 1)
        end
    end

    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)[2:end]
            q[:w, k, i] =
                bound(q[:w, k, i] * aux[:HVSD_w, k, i], params[:w_bounds])
            q[:a, k, i] =
                bound(q[:a, k, i] * aux[:HVSD_w, k, i], params[:a_bounds])

            weight = aux[:HVSD_w, k, i]
            q[:θ_liq, k, i] =
                weight * q[:θ_liq, k, i] + (1 - weight) * q[:θ_liq, k, gm]
            q[:q_tot, k, i] =
                weight * q[:q_tot, k, i] + (1 - weight) * q[:q_tot, k, gm]
        end
    end
end

function compute_cv_gm!(grid, q, ϕ, ψ, cv, tke_factor)
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for k in over_elems(grid)
        Δϕ = q[ϕ, k, en] - q[ϕ, k, gm]
        Δψ = q[ψ, k, en] - q[ψ, k, gm]
        q[cv, k, en] =
            tke_factor * q[:a, k, en] * Δϕ * Δψ + q[:a, k, en] * q[cv, k, en]
        @inbounds for i in ud
            Δϕ = q[ϕ, k, i] - q[ϕ, k, gm]
            Δψ = q[ψ, k, i] - q[ψ, k, gm]
            q[cv, k, en] += tke_factor * q[:a, k, i] * Δϕ * Δψ
        end
    end
end

function compute_mf_gm!(grid, q, aux)
    gm, en, ud, sd, al = allcombinations(q)
    domain_c = over_elems_real(grid)

    @inbounds for i in ud
        @inbounds for k in domain_c
            aux[:mf_aux, k, i] =
                (q[:w, k, i] - q[:w, k, en]) * aux[:ρ_0, k] * q[:a, k, i]
        end
    end

    @inbounds for k in domain_c
        aux[:mf_θ_liq, k] = sum([
            aux[:mf_aux, k, i] * (q[:θ_liq, k, i] - q[:θ_liq, k, en])
            for i in ud
        ])
        aux[:mf_q_tot, k] = sum([
            aux[:mf_aux, k, i] * (q[:q_tot, k, i] - q[:q_tot, k, en])
            for i in ud
        ])
    end

    @inbounds for k in domain_c
        aux[:mf_tend_θ_liq, k] =
            -aux[:α_0, k] * grad(aux[:mf_θ_liq, Cut(k)], grid)
        aux[:mf_tend_q_tot, k] =
            -aux[:α_0, k] * grad(aux[:mf_q_tot, Cut(k)], grid)
    end
end

function compute_cv_shear!(grid::Grid{FT}, q, aux, aux_O2, ϕ, ψ, cv) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    is_tke = cv == :tke
    tke_factor = is_tke ? FT(0.5) : 1
    grad_u = 0
    grad_v = 0
    @inbounds for k in over_elems_real(grid)
        if is_tke
            grad_u = ∇_dn(q[:u, Cut(k), gm], grid)
            grad_v = ∇_dn(q[:v, Cut(k), gm], grid)
            grad_ϕ = ∇_dn(q[ϕ, Cut(k), en], grid)
            grad_ψ = ∇_dn(q[ψ, Cut(k), en], grid)
        else
            grad_ϕ = grad(q[ϕ, Cut(k), en], grid)
            grad_ψ = grad(q[ψ, Cut(k), en], grid)
        end
        ρaK = aux[:ρ_0, k] * q[:a, k, en] * aux[:K_h, k, gm]
        aux_O2[cv][:shear, k] =
            tke_factor * 2 * ρaK * (grad_ϕ * grad_ψ + grad_u^2 + grad_v^2)
    end
end

function compute_cv_interdomain_src!(
    grid::Grid{FT},
    q,
    aux,
    aux_O2,
    ϕ,
    ψ,
    cv,
    tke_factor,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for k in over_elems(grid)
        aux_O2[cv][:interdomain, k] = FT(0)
        @inbounds for i in ud
            Δϕ = q[ϕ, k, i] - q[ϕ, k, en]
            Δψ = q[ψ, k, i] - q[ψ, k, en]
            aux_O2[cv][:interdomain, k] +=
                tke_factor * q[:a, k, i] * (1 - q[:a, k, i]) * Δϕ * Δψ
        end
    end
end

function compute_cv_env!(
    grid::Grid{FT},
    q,
    aux,
    aux_O2,
    ϕ,
    ψ,
    cv,
    tke_factor,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for k in over_elems(grid)
        if q[:a, k, en] > 0
            @inbounds for i in sd
                Δϕ = q[ϕ, k, i] - q[ϕ, k, gm]
                Δψ = q[ψ, k, i] - q[ψ, k, gm]
                q[cv, k, en] -= tke_factor * q[:a, k, i] * Δϕ * Δψ
            end
            q[cv, k, en] = q[cv, k, en] / q[:a, k, en]
        else
            q[cv, k, en] = FT(0)
        end
    end
end

function cleanup_covariance!(grid::Grid{FT}, q) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @inbounds for k in over_elems_real(grid)
        q[:tke, k, en] = q[:tke, k, en] < eps(FT) ? FT(0) : q[:tke, k, en]
    end
end


function construct_tridiag_diffusion_O1!(
    grid::Grid{FT},
    q,
    aux,
    Δt,
    tri_diag,
) where {FT}
    gm, en, ud, sd, al = allcombinations(aux)
    k_1 = first_interior(grid, Zmin())
    k_2 = first_interior(grid, Zmax())
    dzi = grid.Δzi
    @inbounds for k in over_elems_real(grid)
        ρaK_cut = tri_diag[:ρaK, Cut(k)]
        X = aux[:ρ_0, k] * q[:a, k, en] / Δt
        Z = ρaK_cut[1] * dzi * dzi
        Y = ρaK_cut[2] * dzi * dzi
        if k == k_1
            Z = FT(0)
        elseif k == k_2
            Y = FT(0)
        end
        tri_diag[:a, k] = -Z / X
        tri_diag[:b, k] = 1 + Y / X + Z / X
        tri_diag[:c, k] = -Y / X
    end
end

function construct_tridiag_diffusion_O2!(
    grid::Grid{FT},
    q,
    aux,
    params,
    tri_diag,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    Δzi = grid.Δzi
    Δzi2 = grid.Δzi2
    dti = 1 / params[:Δt][1]
    k_1 = first_interior(grid, Zmin())
    k_2 = first_interior(grid, Zmax())

    @inbounds for k in over_elems_real(grid)
        ρ_0_cut = aux[:ρ_0, Cut(k)]
        ae_cut = q[:a, Cut(k), en]
        w_cut = q[:w, Cut(k), en]
        ρa_K = q[:a, Cut(k), en] .* aux[:K_h, Cut(k), gm] .* aux[:ρ_0, Cut(k)]

        D_env = sum([
            ρ_0_cut[2] * q[:a, k, i] * q[:w, k, i] * aux[:ε_model, k, i]
            for i in ud
        ])

        l_mix = max(aux[:l_mix, k, gm], FT(1))
        tke_env = max(q[:tke, k, en], FT(0))

        tri_diag[:a, k] = (-ρa_K[2] * Δzi2)
        tri_diag[:b, k] = (
            ρ_0_cut[2] * ae_cut[2] * dti -
            ρ_0_cut[2] * ae_cut[2] * w_cut[2] * Δzi +
            ρa_K[3] * Δzi2 +
            ρa_K[2] * Δzi2 +
            D_env +
            ρ_0_cut[2] * ae_cut[2] * params[:tke_diss_coeff] * sqrt(tke_env) / l_mix
        )
        tri_diag[:c, k] =
            (ρ_0_cut[3] * ae_cut[3] * w_cut[3] * Δzi - ρa_K[3] * Δzi2)
    end

    tri_diag[:a, k_1] = FT(0)
    tri_diag[:b, k_1] = FT(1)
    tri_diag[:c, k_1] = FT(0)

    tri_diag[:b, k_2] += tri_diag[:c, k_2]
    tri_diag[:c, k_2] = FT(0)
end
