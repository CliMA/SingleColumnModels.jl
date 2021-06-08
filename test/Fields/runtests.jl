using SingleColumnModels: CentField
using ClimateMachineCore.DataLayouts: VF
using ClimateMachineCore: Spaces, Fields

n_elems_real = 10 # number of elements
z_min = 0.0
z_max = 100.0
space = Spaces.warp_mesh(
    Spaces.FaceFiniteDifferenceSpace(z_min, z_max, n_elems_real)
);
FT = Float64
n_up = 2

vals = (
    ρ_gm = FT(0),
    w_gm = FT(0),
    θ_liq_gm = FT(0),
    q_tot_gm = FT(0),
    ρ_en = FT(0),
    w_en = FT(0),
    θ_liq_en = FT(0),
    q_tot_en = FT(0),
    up = ntuple(n_up) do i
        (a = FT(0),
        w = FT(0),
        θ_liq = FT(0),
        q_tot = FT(0))
    end,
)

cent_field = CentField(space, vals)

w_up_i = parent(Fields.field_values(getindex.(cent_field.up, 1)).w)[2]
w_up_i = getindex.(cent_field.up, 1).w[k]
w_up_i = f.up[i].w[k]

nothing
