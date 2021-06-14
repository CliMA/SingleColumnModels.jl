module TriDiagSolvers

include("tridiag_solver_funcs.jl")

export solve_tdma!
export solve_tridiag_wrapper!

using LinearAlgebra
using ..FiniteDifferenceGrids
using ..StateVecs

struct DifussionImplicit end
function TridiagonalMatrix(
    q::StateVec,
    aux::StateVec,
    Δt::T,
    ρ::S,
    a::S,
    K::S,
    grid::Grid{T},
    ::DifussionImplicit,
    i::Int,
) where {S <: Symbol, T}
    #! format: off
    a_ = [   -Δt * grid.Δzi2 *      (aux[ρ, Dual(k)] .* q[a, Dual(k), i] .* aux[K, Dual(k), i])[2] for k in over_elems_real(grid)[1:(end - 1)]]
    b_ = [1 + Δt * grid.Δzi2 * (sum((aux[ρ, Dual(k)] .* q[a, Dual(k), i] .* aux[K, Dual(k), i]))) for k in over_elems_real(grid)]
    c_ = [   -Δt * grid.Δzi2 *      (aux[ρ, Dual(k)] .* q[a, Dual(k), i] .* aux[K, Dual(k), i])[1] for k in over_elems_real(grid)[2:end]]
    #! format: on
    return Tridiagonal(a_, b_, c_)
end

"""
    solve_tdma!(q::StateVec, tendencies::StateVec, aux::StateVec, x::S, ρ::S, a_τp1::S, a_τ::S, K::S, grid::Grid{T}, Δt::T, i::Int=1)

Solve tri-diagonal system using Diffusion-implicit time-marching.
"""
function solve_tdma!(
    q::StateVec,
    tendencies::StateVec,
    aux::StateVec,
    x::S,
    ρ::S,
    a_τp1::S,
    a_τ::S,
    K::S,
    grid::Grid{T},
    Δt::T,
    i::Int = 1,
) where {S <: Symbol, T}
    B = [
        aux[a_τ, k, i] / q[a_τp1, k, i] * q[x, k, i] +
        Δt * tendencies[x, k, i] / (aux[ρ, k] * q[a_τp1, k, i])
        for k in over_elems_real(grid)
    ]
    A = TridiagonalMatrix(q, aux, Δt, ρ, a_τp1, K, grid, DifussionImplicit(), i)
    X = inv(A) * B
    assign_real!(q, x, grid, X, i)
end

"""
    solve_tridiag_wrapper!(grid::Grid, sv::StateVec, ϕ::Symbol, i::Int, tri_diag::StateVec)

Solve the tri-diagonal system, given by state vector `tri_diag` for variable `ϕ` in sub-domain `i`.
"""
function solve_tridiag_wrapper!(
    grid::Grid,
    sv::StateVec,
    ϕ::Symbol,
    i::Int,
    tri_diag::StateVec,
)
    f = [tri_diag[:f, k] for k in over_elems_real(grid)]
    a = [tri_diag[:a, k] for k in over_elems_real(grid)[2:end]]
    b = [tri_diag[:b, k] for k in over_elems_real(grid)]
    c = [tri_diag[:c, k] for k in over_elems_real(grid)[1:end-1]]
    A = Tridiagonal(a, b, c)
    x = inv(A) * f
    assign_real!(sv, ϕ, grid, x, i)
end

end
