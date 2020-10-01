using Test

using LinearAlgebra

using SingleColumnModels.TriDiagSolvers
TDMA = TriDiagSolvers

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
    # @show eps(Float64)
    # @show "Full", err
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

