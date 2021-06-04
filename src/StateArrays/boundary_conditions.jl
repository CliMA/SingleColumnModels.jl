# module BoundaryConditions
# using ..VariableTemplates
# using ..BalanceLaws: Prognostic, vars_state

abstract type AbstractBCType end
struct Dirichlet <: AbstractBCType end
struct Neumann <: AbstractBCType end
struct Robin <: AbstractBCType end

abstract type AbstractBCForm{BCT<:AbstractBCType} end; export AbstractBCForm
struct NonHomogeneousBC{BCT,VT,FT} <: AbstractBCForm{BCT}; val::VT; end; export NonHomogeneousBC
struct HomogeneousBC{BCT,FT<:AbstractFloat} <: AbstractBCForm{BCT} end; export HomogeneousBC

bc_given_form(bc, ::Type{HomogeneousBC}) = bc.h
bc_given_form(bc, ::Type{NonHomogeneousBC}) = bc.nh

export BCValType
struct BCValType{V} end
get_value(::BCValType{Symbol}, ::ZBoundary) = MethodError("No BC found")

export BCs
struct BCs{HT,NHT,L}
    h::HT
    nh::NHT
    loc::L
    function BCs(::Type{FT},
            var::Symbol;
            loc=Zmin(),
            bctype=Dirichlet(),
            val=BCValType{var}()
        ) where {FT}
        T = typeof(bctype)
        VT = typeof(val)
        L = typeof(loc)
        h = HomogeneousBC{T,FT}()
        nh = NonHomogeneousBC{T,VT,FT}(val)
        return new{typeof(h),typeof(nh),L}(h, nh)
    end
end

val(bc::HomogeneousBC{BCT,FT}, zb) where {BCT,FT} = FT(0)
val(bc::NonHomogeneousBC{BCT,VT,FT}, zb) where {BCT,VT,FT} = FT(get_value(bc.val, zb))

export apply_bcs!
function apply_bcs!(data, grid, bcs::NamedTuple,
    bcform::Union{Type{HomogeneousBC},Type{NonHomogeneousBC}})
    for bc_generic in bcs
        bc = bc_given_form(bc_generic, bcform)
        set_bcs!(data, grid, bc, bc_generic.loc)
    end
end

# function apply_bcs!(data, grid,
#     bcform::Union{Type{HomogeneousBC},Type{NonHomogeneousBC}})
#     # vs = vars_state(state.model, Prognostic(), eltype(state))
#     # bcs = state.model.bcs
#     # FT = eltype(grid)
#     # for ftc in flattened_tup_chain(vs)
#     #     sym = ftc[end]
#     #     vi = varsindex(vs, ftc...)
#     #     progstate = view(state.data, :, vi)
#     #     if hasproperty(bcs, sym)
#     #         apply_bcs_single_var!(progstate, grid, bc_given_form(getfield(bcs, sym), bcform, Zmin()), Zmin())
#     #         apply_bcs_single_var!(progstate, grid, bc_given_form(getfield(bcs, sym), bcform, Zmax()), Zmax())
#     #     end
#     # end
#     FT = eltype(grid)
#     loc = bc.loc
#     set_bcs!(data, grid, bc_given_form(bc, bcform), loc)
# end
set_bcs!(prog, grid, bc::AbstractBCForm{BCT}, zb::Zmin) where {BCT<:Dirichlet} = (prog[1] = 2*val(bc, zb)-prog[2])
set_bcs!(prog, grid, bc::AbstractBCForm{BCT}, zb::Zmax) where {BCT<:Dirichlet} = (prog[end] = 2*val(bc, zb)-prog[end-1])
set_bcs!(prog, grid, bc::AbstractBCForm{BCT}, zb::Zmin) where {BCT<:Neumann} = (prog[1] = prog[2] + n_hat(grid, zb)*val(bc, zb)*grid.Δz)
set_bcs!(prog, grid, bc::AbstractBCForm{BCT}, zb::Zmax) where {BCT<:Neumann} = (prog[end] = prog[end-1] + n_hat(grid, zb)*val(bc, zb)*grid.Δz)

# end