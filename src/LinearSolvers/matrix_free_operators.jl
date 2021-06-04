#### Matrix-free operators

export MatrixFree∇²
struct MatrixFree∇²{FT,G,ABC,BCS,AT1} <: AbstractMatrix{FT}
    grid::G
    apply_bcs!::ABC
    bcs::BCS
    vol::AT1
end
function MatrixFree∇²(grid,apply_bcs!,bcs)
    FT = eltype(grid)
    vol = similar(grid.zc)
    vol .= grid.Δz
    args = (grid,apply_bcs!,bcs,vol)
    return MatrixFree∇²{FT,typeof.(args)...}(args...)
end

export compute_Ax_bc!
function compute_Ax_bc!(pcg, grid, x)
    @unpack_fields pcg A Ax_bc x_bc
    @unpack_fields A apply_bcs!
    x_bc .= 0
    apply_bcs!(x_bc, grid, A.bcs.x, NonHomogeneousBC)
    mul!(Ax_bc, A, x_bc;bc_form=NonHomogeneousBC)
    for i in over_elems_ghost(grid) # assuming not periodic
        Ax_bc[i] = x[i]
    end
    return nothing
end

"""
Laplacian operator:
    `` vol*A = vol*∇•∇``
    `` A = ∇•∇`` if `scale=false`
"""
function LinearAlgebra.mul!(
    Ax::AbstractVector,A::MatrixFree∇²,
    x::AbstractVector;scale=true,bc_form=HomogeneousBC)
    @unpack_fields A grid vol bcs
    A.apply_bcs!(x, grid, bcs.x, bc_form)
    @inbounds for i in over_elems_real(grid)
        Ax[i] = (x[i-1]+x[i+1]-2*x[i])/grid.Δz
    end
    scale && (Ax .*= vol)
    return nothing
end
