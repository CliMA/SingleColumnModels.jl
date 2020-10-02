module ConjugateGradientMethods

using ..Utilities
using ..FiniteDifferenceGrids
using LinearAlgebra

export ∇², solve!

include("cg_types.jl")

all_Neumann(::PCG) = false
function subtract_physical_mean! end

"""
    ∇²(∇²x,x,grid,params)
Laplacian operator: ``A = ∇²•(∇)x``
"""
function ∇²(∇²x,x,grid,params)
    @inbounds for i in over_elems_real(grid)
        ∇²x[i] = (x[i-1]+x[i+1]-2*x[i])/grid.Δz
    end
end

function modify_RHS!(pcg, grid, x, b) end

Base.@kwdef mutable struct Norms{FT}
    L₂::FT = FT(0)
    L₁::FT = FT(0)
    L∞::FT = FT(0)
end
norms(x) = (L₂ = norm(x), L₁ = norm(x, 1), L∞ = norm(x, Inf))

function solve!(x, b, pcg, grid)
    modify_RHS!(pcg, grid, x, b)
    @unpack_fields pcg operator! apply_BCs! isp mfp Ax k z p r M⁻¹ vol tempk
    # ********************* START PCG ALGORITHM *********************
    res_norm₀ = (L₂ = norm(r), L₁ = norm(r, 1), L∞ = norm(r, Inf))
    isp.iter_performed = 0
    isnan(res_norm₀.L₂) && error("res_norm₀.L₂ is NaNs")
    z = M⁻¹ .* r
    p = z
    ρ_k = r'*z
    res_norm = (L₂=sqrt(abs(ρ_k)),)
    if !exit_loop(isp) # Only do PCG if necessary!
        for i=1:iter_max(isp)
            apply_BCs!(p, grid)
            operator!(Ax,p,grid,(mfp=mfp,k=k,tempk=tempk))
            Ax .*= vol
            α = ρ_k/(p'*Ax)
            x .= x .+ α * p
            apply_BCs!(x, grid) # Needed for PPE
            update_iter!(isp)
            r .= r .- α*Ax # r = r - α Ap
            # if check_res(isp)
                res_norm = (L₂=sqrt(abs(r'*r)),)
                @show res_norm
                update_exit_loop!(isp,res_norm.L₂,res_norm₀.L₂)
                exit_loop(isp) && break
            # end
            z .= M⁻¹ .* r
            ρ_kp1 = z'*r
            p .= p * ρ_kp1/ρ_k .+ z # p = p β + z
            ρ_k = ρ_kp1
            update_exit_loop!(isp)
        end
    else
      apply_BCs!(x, grid) # BCs should still be enforced!
      update_iter!(isp)
    end
    @info res_norm
    isnan(res_norm.L₂) && error("res_norm is NaNs (after solve)")
    return nothing
end

end