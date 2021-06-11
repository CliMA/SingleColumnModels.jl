###### Pressure models

abstract type PressureModel end

Base.@kwdef struct SCAMPyPressure{FT} <: PressureModel
    buoy_coeff::FT = FT(0)
    drag_coeff::FT = FT(0)
    plume_spacing::FT = FT(0)
end

Base.@kwdef struct NormalMode{FT} <: PressureModel
    nh_avd_coef::FT = FT(0)
    buoy_coef::FT = FT(0)
    nh_drag_coef::FT = FT(0)
end

"""
    compute_pressure!

Define pressure field

 - `aux[:nh_press, k, i]`

for all `k` and all `i`
"""
function compute_pressure! end

function compute_pressure!(
    grid::Grid{FT},
    q,
    aux,
    params,
    model::SCAMPyPressure,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    p_coeff = model.drag_coeff / model.plume_spacing

    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            a_k = q[:a, k, i]
            w_env = q[:w, k, en]
            ρ_k = aux[:ρ_0, k]
            w_i = q[:w, k, i]
            B_k = aux[:buoy, k, i]
            ρa_k = ρ_k * a_k

            press_buoy = -ρa_k * B_k * model.buoy_coeff
            press_drag = -ρa_k * (p_coeff * (w_i - w_env)^2 / sqrt(a_k))
            aux[:nh_press, k, i] = press_buoy + press_drag
        end
    end
end

function compute_pressure!(
    grid::Grid{FT},
    q,
    aux,
    params,
    model::NormalMode,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    @unpack UpdVar = params

    @inbounds for i in ud
        @inbounds for k in over_elems_real(grid)
            a_k = q[:a, k, i]
            w_env = q[:w, k, en]
            ρ_k = aux[:ρ_0, k]
            w_i = q[:w, k, i]
            # dwdz = ∇_z_upwind(q[:w, k-1:k+1, i], q[:w, k-1:k+1, i], grid)
            dwdz = (q[:w, k + 1, i] - q[:w, k - 1, i]) * grid.Δzi * 0.5
            B_k = aux[:buoy, k, i]
            ρa_k = ρ_k * a_k

            nh_press_buoy = -ρa_k * B_k * model.buoy_coef
            nh_pressure_adv = ρa_k * model.nh_avd_coef * w_i * dwdz
            nh_pressure_drag =
                -ρa_k * model.nh_drag_coef * (w_i - w_env) * abs(w_i - w_env) /
                max(UpdVar[i].cloud.updraft_top, 500.0)

            aux[:nh_press, k, i] =
                nh_press_buoy + nh_pressure_adv + nh_pressure_drag
        end
    end
end

function compute_tke_pressure!(
    grid::Grid{FT},
    q,
    aux,
    aux_O2,
    cv,
    params,
) where {FT}
    gm, en, ud, sd, al = allcombinations(q)
    is_tke = cv == :tke
    @inbounds for k in over_elems_real(grid)
        aux_O2[cv][:press, k] = FT(0)
        @inbounds for i in ud
            aux_O2[cv][:press, k] = FT(0)
            if is_tke
                @inbounds for i in ud
                    w_env = q[:w, k, en]
                    w_i = q[:w, k, i]
                    aux_O2[cv][:press, k] +=
                        (w_env - w_i) * aux[:nh_press, k, i]
                end
            end
        end
    end
end
