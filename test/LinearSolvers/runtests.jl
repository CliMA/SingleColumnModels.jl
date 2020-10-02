using Test

using LinearAlgebra

using SingleColumnModels.TriDiagSolvers
TDMA = TriDiagSolvers

# @testset "TriDiagSolvers" begin
#   N = 2:10
#   for n in N
#     dl = rand(n-1)
#     du = rand(n-1)
#     d = rand(n);

#     A = Array(Tridiagonal(dl, d, du))
#     b = rand(length(d))

#     x_correct = inv(A)*b

#     xtemp = zeros(n)
#     gamma = zeros(n-1)
#     beta = zeros(n)
#     x_TDMA = zeros(n)

#     solve_tridiag!(x_TDMA, b, dl, d, du, n, xtemp, gamma, beta)
#     tol = eps(Float32)

#     err = [abs(x-y) for (x, y) in zip(x_correct, x_TDMA)]
#     # @show eps(Float64)
#     # @show "Full", err
#     if !all([x<tol for x in err])
#         @show tol
#         @show "Full", err
#     end
#     @test all([x<tol for x in err])

#     init_β_γ!(beta, gamma, dl, d, du, n)
#     solve_tridiag_stored!(x_TDMA, b, dl, beta, gamma, n, xtemp)

#     err = [abs(x-y) for (x, y) in zip(x_correct, x_TDMA)]
#     # @show "Stored", err
#     @test all([x<tol for x in err])
#   end
# end

using SingleColumnModels.ConjugateGradientMethods
using SingleColumnModels.Utilities
import SingleColumnModels.ConjugateGradientMethods: modify_RHS!
using SingleColumnModels.FiniteDifferenceGrids

struct Implicit end
struct Explicit end
function apply_BCs!(x, grid)
    x[1] = -x[2]
    x[end] = -x[end-1]
end
function apply_BCs!(x, grid, implicit) # (i.e., homogeneous BCs)
    x[1] = -x[2]
    x[end] = -x[end-1]
end
all_Neumann(::PCG) = false
function modify_RHS!(pcg, grid, x, b)
    @unpack_fields pcg operator! apply_BCs! mfp Ax x_bc k r vol tempx tempk
    r .= b
    Ax_bc = similar(x)
    all_Neumann(pcg) && subtract_physical_mean!(r,vol) # Not sure correct loc
    # compute x_bc/Ax_bc:
    x_bc .= 0
    apply_BCs!(x_bc, grid)
    operator!(Ax_bc,x_bc,grid,(mfp=mfp,k=k,tempk=tempk))
    for i in over_elems_ghost(grid) # assuming not periodic
        Ax_bc[i] = x[i]
        vol[i] = 0
    end
    r .= r .- Ax_bc                  # r = (b_mod - Ax_bc - Ax)
    apply_BCs!(x, grid, Implicit()) # to make operator implicit
    operator!(Ax,x,grid,(mfp=mfp,k=k,tempk=tempk))
    r .= (r .- Ax) .* vol        # r = vol*(b_mod - Ax_bc - Ax_mod)
end
using DelimitedFiles
@testset "CG" begin
    FT = Float64
    mfp = MatrixFreeParams{FT}()
    isp = IterativeSolverParams{FT}(;
        exit_cond = ExitCondition{FT}(;iter_max=20))
    grid = UniformGrid(0.0, 1.0, 100)
    p = 2FT(π)
    x = sin.(p * grid.zc)
    x_exact = similar(x)
    x_exact .= x

    prec!(M⁻¹,grid,k,coeff,tempx) = (M⁻¹ .= 1)
    pcg = PCG(∇², prec!, apply_BCs!, grid, x;isp=isp)
    b = similar(x)
    b .= 0
    ∇²(b,x,grid,(;))
    # b .= -p^2*sin.(p * grid.zc)
    @show b
    x .= 0
    solve!(x, b, pcg, grid)
    @show b
    @show pcg.isp
    zn = grid.ze[2:end-1]
    apply_BCs!(x, grid)
    @show x
    xn = (x[1:end-1]+x[2:end])/2
    xn_exact = (x_exact[1:end-1]+x_exact[2:end])/2
    bn = (b[1:end-1]+b[2:end])/2
    writedlm("z.csv", zn, ',')
    writedlm("x.csv", xn, ',')
    writedlm("x_exact.csv", xn_exact, ',')
    writedlm("b.csv", bn, ',')
    err = norm(abs.(x .- x_exact))
    @test err < sqrt(eps(FT))/10
    @show sqrt(eps(FT))/10
end

