
export IterativeSolverParams,
    ExitCondition,
    TimeSplit,
    MatrixFreeParams,
    PCG

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

struct PCG{OP,P,ABC,AT1,AT2,ISP,MFP}
    operator!::OP
    prec!::P
    apply_BCs!::ABC
    tempx::AT1
    p::AT1
    x_bc::AT1
    r::AT1
    Ax::AT1
    vol::AT1
    M⁻¹::AT1
    z::AT1
    k::AT2
    tempk::AT2
    isp::ISP
    mfp::MFP
    function PCG(
            operator!,
            prec!,
            apply_BCs!,
            grid,
            x::AbstractArray{FT};
            isp = IterativeSolverParams{FT}(),
            mfp = MatrixFreeParams{FT}(),
            k = similar(x)) where {FT}
        p = similar(x) # Copies BCs for x
        x_bc = similar(x) # Copies BCs for x
        tempx = similar(x)
        r = similar(tempx)
        Ax = similar(tempx)
        vol = similar(tempx)
        z = similar(tempx)
        M⁻¹ = similar(tempx)
        k = similar(k)
        tempk = similar(k)
        vol .= grid.Δz
        prec!(M⁻¹,grid,k,mfp.coeff.implicit,tempx) # mfp%coeff_implicit reasonable estimate
        args = (operator!,prec!,apply_BCs!,x,k,isp,mfp)
        OP,P,ABC,AT1,AT2,ISP,MFP = typeof.(args)
        return new{OP,P,ABC,AT1,AT2,ISP,MFP}(
            operator!,prec!,apply_BCs!,tempx,p,
            x_bc,r,Ax,vol,M⁻¹,z,k,tempk,isp,mfp)
    end
end

