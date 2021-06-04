module ConjugateGradientMethods

using ..Utilities
using ..FiniteDifferenceGrids
using ..StateArrays
using LinearAlgebra

export solve!

include("cg_types.jl")
include("matrix_free_operators.jl")

function solve!(x, b, solver::PCG, grid)
    @unpack_fields solver A isp Ax z p r M⁻¹
    r .= b
    res_norm₀ = (L₂ = norm(r), L₁ = norm(r, 1), L∞ = norm(r, Inf))
    isp.iter_performed = 0
    isnan(res_norm₀.L₂) && error("res_norm₀.L₂ is NaNs")
    z = M⁻¹ .* r
    p = z
    ρ_k = r'*z
    res_norm = (L₂=sqrt(abs(ρ_k)),)
    if !exit_loop(isp) # Only iterate if necessary!
        for i=1:iter_max(isp)
            mul!(Ax,A,p)
            α = ρ_k/(p'*Ax)
            x .+= α * p
            update_iter!(isp)
            r .-= α*Ax
            if check_res(isp)
                res_norm = (L₂=sqrt(abs(r'*r)),)
                @show res_norm
                update_exit_loop!(isp,res_norm.L₂,res_norm₀.L₂)
                exit_loop(isp) && break
            end
            z .= M⁻¹ .* r
            ρ_kp1 = z'*r
            p .= p * ρ_kp1/ρ_k .+ z # p = p β + z
            ρ_k = ρ_kp1
            update_exit_loop!(isp)
        end
    end
    @info res_norm
    isnan(res_norm.L₂) && error("res_norm is NaNs (after solve)")
    return nothing
end

end