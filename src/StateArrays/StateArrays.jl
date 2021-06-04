module StateArrays

using ..FiniteDifferenceGrids
# using ..BalanceLaws: Prognostic, AbstractStateType

export State
struct State{FT,N,AFT,M} <: AbstractArray{FT,N}
    model::M
    data::AFT
end

# function State(model, grid)
#     vs = vars_state(model, Prognostic(), eltype(grid))
#     N = varsize(vs)
#     FT = eltype(grid)
#     data = zeros(FT, grid.n_elem, N)
#     return State{FT, 2, typeof(data), typeof(model)}(model, data)
# end

# Base.size(s::State) = Base.size(s.data)
# Base.getindex(s::State, i) = Base.getindex(s.data, i)
# Base.getindex(s::State, i, j) = Base.getindex(s.data, i, j)
# Base.setindex!(s::State, i, v) = Base.setindex!(s.data, i, v)
# Base.setindex!(s::State, i, j, v) = Base.setindex!(s.data, i, j, v)
# get_N(state::State{FT,N}) where {FT,N} = N
# Base.similar(state::State) =
#     State{eltype(state),
#           get_N(state),
#           typeof(state.data),
#           typeof(state.model)}(
#       state.model,
#       similar(state.data)
#     )

# export field
# field(state::State, st::AbstractStateType, syms...) =
#     view(state.data, :, varsindex(vars_state(state.model, st, eltype(state)), syms...))
# field(state::State, sym, st::AbstractStateType = Prognostic()) =
#     field(state, st, sym)

# # Base.getproperty(state::State, syms::Symbol) =
# #     field(state, Prognostic(), syms)

# # Base.getproperty(state::State, syms::Tuple) =
# #     field(state, Prognostic(), syms...)

# export PrognosticFieldType
# Base.@kwdef struct PrognosticFieldType{BC}
#     bcs::BC
# end
# export DiagnosticFieldType
# Base.@kwdef struct DiagnosticFieldType{S} end

include("boundary_conditions.jl")

end