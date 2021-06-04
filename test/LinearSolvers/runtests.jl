using Test
using LinearAlgebra
using SingleColumnModels.TriDiagSolvers
const TDMA = TriDiagSolvers

@testset "TriDiagSolvers" begin
  N = 2:10
  for n in N
    dl = rand(n-1)
    du = rand(n-1)
    d = rand(n);

    A = Array(Tridiagonal(dl, d, du))
    b = rand(length(d))

    x_correct = inv(A)*b

    xtemp = zeros(n)
    gamma = zeros(n-1)
    beta = zeros(n)
    x_TDMA = zeros(n)

    solve_tridiag!(x_TDMA, b, dl, d, du, n, xtemp, gamma, beta)
    tol = eps(Float32)

    err = [abs(x-y) for (x, y) in zip(x_correct, x_TDMA)]

    if !all([x<tol for x in err])
        @show tol
        @show "Full", err
    end
    @test all([x<tol for x in err])

    init_β_γ!(beta, gamma, dl, d, du, n)
    solve_tridiag_stored!(x_TDMA, b, dl, beta, gamma, n, xtemp)

    err = [abs(x-y) for (x, y) in zip(x_correct, x_TDMA)]
    # @show "Stored", err
    @test all([x<tol for x in err])
  end
end

using DelimitedFiles
using SingleColumnModels.ConjugateGradientMethods
using SingleColumnModels.FiniteDifferenceGrids
using SingleColumnModels.Utilities
using SingleColumnModels.StateArrays
const SA = StateArrays
import SingleColumnModels.StateArrays: get_value

get_value(::BCValType{:x}, ::Zmin) = 0
get_value(::BCValType{:x}, ::Zmax) = 0
all_Neumann(::PCG) = false
function subtract_physical_mean! end
function modify_RHS!(pcg, grid, x, b)
    @unpack_fields pcg A Ax_bc Ax r
    @unpack_fields A apply_bcs!
    compute_Ax_bc!(pcg, grid, x)
    r .= b .- Ax_bc       # r = (b_mod - Ax_bc - Ax)
    mul!(Ax, A, x)
    r .= (r .- Ax)        # r = (b_mod - Ax_bc - Ax_mod)
    b .= r
end

@testset "∇²u = f = ∇²(sin(2πx)) = -2π²sin(2πx), u = 0 ∈ ∂Ω" begin
    FT = Float64
    isp = IterativeSolverParams{FT}(;
        exit_cond = ExitCondition{FT}(;iter_max=20))
    grid = UniformGrid(0.0, 1.0, 100)
    p = 2FT(π)
    x = sin.(p * grid.zc)
    x_exact = similar(x)
    x_exact .= x
    bcs = (
           x=(zmin=BCs(FT,:x;loc=Zmin(),bctype=SA.Dirichlet()),
              zmax=BCs(FT,:x;loc=Zmax(),bctype=SA.Dirichlet())
              ),
          )

    A = MatrixFree∇²(grid, apply_bcs!, bcs)
    pcg = PCG(A, grid, x; isp=isp)
    b = similar(x)
    b .= 0
    mul!(b,A,x)
    x .= 0

    modify_RHS!(pcg, grid, x, b)
    solve!(x, b, pcg, grid)
    # BCs should still be enforced!
    apply_bcs!(x, grid, A.bcs.x, NonHomogeneousBC)

    zn = grid.ze[2:end-1]
    xn = (x[1:end-1]+x[2:end])/2
    xn_exact = (x_exact[1:end-1]+x_exact[2:end])/2
    bn = (b[1:end-1]+b[2:end])/2
    writedlm("z1.csv", zn, ',')
    writedlm("x1.csv", xn, ',')
    writedlm("x1_exact.csv", xn_exact, ',')
    writedlm("b1.csv", bn, ',')
    err = norm(abs.(x .- x_exact))
    @test err < sqrt(eps(FT))/10
    @show sqrt(eps(FT))/10
end

get_value(::BCValType{:x}, ::Zmin) = 1
get_value(::BCValType{:x}, ::Zmax) = 0
@testset "∇²u = f, u = (0, 1)" begin
    FT = Float64
    isp = IterativeSolverParams{FT}(;
        exit_cond = ExitCondition{FT}(;iter_max=100))
    grid = UniformGrid(0.0, 1.0, 100)
    p = 2FT(π)
    x = sin.(p * grid.zc)
    x_exact = similar(x)
    x_exact .= 1 .- grid.zc
    bcs = (
           x=(zmin=BCs(FT,:x;loc=Zmin(),bctype=SA.Dirichlet()),
              zmax=BCs(FT,:x;loc=Zmax(),bctype=SA.Dirichlet())
              ),
          )

    A = MatrixFree∇²(grid, apply_bcs!, bcs)
    pcg = PCG(A, grid, x; prec! = diagonal_preconditioner(grid), isp=isp)
    b = similar(x)
    b .= 0
    x .= 0

    modify_RHS!(pcg, grid, x, b)
    solve!(x, b, pcg, grid)
    # BCs should still be enforced!
    apply_bcs!(x, grid, A.bcs.x, NonHomogeneousBC)

    zn = grid.ze[2:end-1]
    xn = (x[1:end-1]+x[2:end])/2
    xn_exact = (x_exact[1:end-1]+x_exact[2:end])/2
    bn = (b[1:end-1]+b[2:end])/2
    writedlm("z2.csv", zn, ',')
    writedlm("x2.csv", xn, ',')
    writedlm("x2_exact.csv", xn_exact, ',')
    writedlm("b2.csv", bn, ',')
    if isp.exit_cond.iter_max == 100
        x_iter_100_diag_prec = [1.0147629510200533, 0.9852370489799467, 0.9569376782527182, 0.9276286065754576, 0.898685970896808, 0.8698520532101993, 0.8400751396459679, 0.8117410520481221, 0.7814859999122173, 0.7532205341162982, 0.7233427485381876, 0.6946820523228039, 0.6658226000993187, 0.6368446315893351, 0.6092211421803455, 0.5805419614032142, 0.5539197844837326, 0.5263651625403525, 0.5005800045179679, 0.47461353490026537, 0.44975831853681525, 0.42566651148128526, 0.4019735609443584, 0.3797499432037947, 0.3575840491253121, 0.3370503286130786, 0.3166984735229828, 0.2975800013145602, 0.2792418755951479, 0.2614964121737386, 0.24496248239027366, 0.2288199412073564, 0.21377832903150726, 0.1993253679200035, 0.18557084146377842, 0.1726583699642897, 0.16020995653160452, 0.14850833072534272, 0.13733317225316546, 0.1266729111583249, 0.11654964680698653, 0.10681256471162176, 0.09753318965566063, 0.08863984947846325, 0.08010890369088519, 0.07199005054566907, 0.06425483517235935, 0.05692452347295669, 0.05002226303397975, 0.04355710677682626, 0.03755352574037485, 0.03202238368979958, 0.02697934548418765, 0.02243818513192596, 0.018393522650310128, 0.014838414895836307, 0.011761073699617743, 0.00913616026817904, 0.006935054432330808, 0.00512820441596293, 0.003680339826513149, 0.002552214416470799, 0.0017025844735885334, 0.0010874825946037736, 0.0006615113612751936, 0.00038098624554370927, 0.00020648090304070514, 0.00010459857095372124, 4.913734105431334e-5, 2.1214232519249162e-5, 8.340199106826197e-6, 2.961098264640345e-6, 9.430651087255557e-7, 2.68088402548102e-7, 6.778566007827917e-8, 1.5207642817253702e-8, 3.0218985979148564e-9, 5.310842465473636e-10, 8.243547309596993e-11, 1.1284483099147368e-11, 1.359890879619392e-12, 1.4396712248275018e-13, 1.3355200324647745e-14, 1.08226845998989e-15, 7.633728940504137e-17, 4.666537402523848e-18, 2.4599875637176523e-19, 1.1117766018248825e-20, 4.278622273065726e-22, 1.3911279821816724e-23, 3.7862490232508872e-25, 8.533643197664579e-27, 1.5724417130438639e-28, 2.3325708140599975e-30, 2.73342414067111e-32, 2.470780543426337e-34, 1.6694902493499245e-36, 8.069662229400939e-39, 2.6077495818044567e-41, 4.992278584783237e-44, 4.2330028563049355e-47, -4.2330028563049355e-47]
        err = norm(abs.(x .- x_iter_100_diag_prec))
        @test err < sqrt(eps(FT))
    end
end

