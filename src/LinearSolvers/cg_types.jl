
export IterativeSolverParams,
    ExitCondition,
    TimeSplit,
    MatrixFreeParams,
    PCG

Base.@kwdef mutable struct Norms{FT}
    L₂::FT = FT(0)
    L₁::FT = FT(0)
    L∞::FT = FT(0)
end
norms(x) = (L₂ = norm(x), L₁ = norm(x, 1), L∞ = norm(x, Inf))

abstract type AbstractPreconditioner end

struct Diagonal <: AbstractPreconditioner end

struct TimeSplit{FT}
    "coefficient of explicit terms without time discretization"
    explicit::FT
    "coefficient of implicit terms without time discretization"
    implicit::FT
end

Base.@kwdef mutable struct MatrixFreeParams{FT}
    "weight of implicit treatment (1 = Backward Euler, .5 = Crank Nicholson)"
    α::FT = FT(0.5)
    "weight of explicit treatment"
    β::FT = 1 - α
    "coefficient of terms in non-discretized equation"
    coeff_natural::FT = FT(0) # e.g., 1/Reynolds
    "Δt*coeff.implicit/coeff_unsteady"
    coeff_implicit_time_split::FT = FT(0)
    "coefficients without time discretization"
    coeff::TimeSplit{FT} = TimeSplit{FT}(coeff_natural*β, -coeff_natural*α)
end

update_mfp!(mfp, Δt) = (mfp.coeff_implicit_time_split = Δt*mfp.coeff.implicit)

Base.@kwdef mutable struct ExitCondition{FT}
    iter_max::Int = 10
    tol_abs::FT = sqrt(eps(FT))
    tol_rel::FT = sqrt(eps(FT))
end

Base.@kwdef mutable struct IterativeSolverParams{FT}
    iter_performed::Int = 0
    exit_cond = ExitCondition{FT}()
    n_skip_check_res::Int = 5
    exit_loop::Tuple{Bool,Bool,Bool} = (false,false,false)
end
function Base.show(io::IO, isp::IterativeSolverParams)
  println(io, "\n----------------------- IterativeSolverParams")
  println(io, "iter_performed   = $(isp.iter_performed)")
  println(io, "iter_max         = $(isp.exit_cond.iter_max)")
  println(io, "tol_abs          = $(isp.exit_cond.tol_abs)")
  println(io, "tol_rel          = $(isp.exit_cond.tol_rel)")
  println(io, "n_skip_check_res = $(isp.n_skip_check_res)")
  println(io, "exit_loop        = $(isp.exit_loop)")
  println(io, "\n-----------------------")
end
iter_max(isp::IterativeSolverParams) = isp.exit_cond.iter_max
exit_loop(isp::IterativeSolverParams) = any(isp.exit_loop)

function update_iter!(isp::IterativeSolverParams)
    isp.iter_performed += 1
end

function update_exit_loop!(isp,res,res₀)
    @unpack_fields isp exit_cond
    EL_1 = res/res₀ < exit_cond.tol_rel
    EL_2 = res < exit_cond.tol_abs
    EL_3 = isp.iter_performed >= exit_cond.iter_max
    isp.exit_loop = (EL_1,EL_2,EL_3)
end

function update_exit_loop!(isp)
    @unpack_fields isp exit_cond
    isp.exit_loop = (isp.exit_loop[1:2]...,isp.iter_performed>exit_cond.iter_max)
end

check_res(isp::IterativeSolverParams) = mod(isp.iter_performed,isp.n_skip_check_res) == 0

export diagonal_preconditioner
diagonal_preconditioner(grid) = x -> fill!(x, 1/grid.Δz^2)

struct PCG{AT2,P,AT1,VT2,ISP,MFP}
    A::AT2
    prec!::P
    Ax_bc::AT1
    p::AT1
    x_bc::AT1
    r::AT1
    Ax::AT1
    M⁻¹::AT1
    z::AT1
    k::VT2
    tempk::VT2
    isp::ISP
    mfp::MFP
    function PCG(
            A,
            grid,
            x::AbstractArray{FT};
            prec! = x->fill!(x, 1),
            isp = IterativeSolverParams{FT}(),
            mfp = MatrixFreeParams{FT}(),
            k = similar(x)) where {FT}
        p = similar(x) # Copies BCs for x
        x_bc = similar(x) # Copies BCs for x
        Ax_bc = similar(x)
        r = similar(Ax_bc)
        Ax = similar(Ax_bc)
        z = similar(Ax_bc)
        M⁻¹ = similar(Ax_bc)
        k = similar(k)
        tempk = similar(k)
        prec!(M⁻¹) # mfp.coeff_implicit reasonable estimate
        args = (A,prec!,x,k,isp,mfp)
        AT2,P,AT1,VT2,ISP,MFP = typeof.(args)
        return new{AT2,P,AT1,VT2,ISP,MFP}(
            A,prec!,Ax_bc,p,
            x_bc,r,Ax,M⁻¹,z,k,tempk,isp,mfp)
    end
end

