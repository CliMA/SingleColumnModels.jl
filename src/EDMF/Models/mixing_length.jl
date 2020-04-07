##### Mixing length models

abstract type MixingLengthModel end

struct ConstantMixingLength{FT} <: MixingLengthModel
  value::FT
end

function compute_mixing_length!(grid::Grid{FT}, q, tmp, params, model::ConstantMixingLength) where FT
  gm, en, ud, sd, al = allcombinations(q)
  @inbounds for k in over_elems(grid)
    tmp[:l_mix, k, gm] = model.value
  end
end
