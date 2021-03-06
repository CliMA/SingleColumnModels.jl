using Test, Printf

using NCDatasets
using SingleColumnModels.DomainDecompositions
using SingleColumnModels.FiniteDifferenceGrids
using SingleColumnModels.StateVecs

n_elems_real = 10 # number of elements

grid = UniformGrid(0.0, 1.0, n_elems_real)
FT = eltype(grid)
vars = (
    (:ρ_0, DomainSubSet(gm = true)),
    (:a, DomainSubSet(gm = true, en = true, ud = true)),
    (:w, DomainSubSet(gm = true, en = true, ud = true)),
)
domain_set = DomainSet(gm = 1, en = 1, ud = 2)
q = StateVec(vars, grid, domain_set)
q_compare = deepcopy(q)
gm, en, ud, sd, al = allcombinations(q)

@testset "State Vectors" begin
    for k in over_elems(grid)
        q[:ρ_0, k] = 2.0
    end
    @test all([q[:ρ_0, k] for k in over_elems(grid)] .== 2.0)

    for k in over_elems(grid)
        q[:w, k, gm] = 3.0
    end
    @test all([q[:w, k, gm] for k in over_elems(grid)] .== 3.0)

    @test_throws BoundsError q[:w, 1, 1000] = 3.0
    @test_throws BoundsError q[:ρ_0, 1, en] = 3.0

    for k in over_elems(grid)
        ρ_0_e = q[:ρ_0, k]
        for i in alldomains(q)
            w_0_e_i = q[:w, k, i]
        end
    end
    sprint(show, q)

    assign!(q, grid, 0.0)
    k = 1
    apply_Dirichlet!(q, :ρ_0, grid, FT(2), Zmin())
    @test q[:ρ_0, k] ≈ 4
    apply_Neumann!(q, :ρ_0, grid, -2 / grid.Δz, Zmin())
    @test q[:ρ_0, k] ≈ 2

    k = grid.n_elem
    apply_Dirichlet!(q, :ρ_0, grid, FT(2), Zmax())
    @test q[:ρ_0, k] ≈ 4
    apply_Neumann!(q, :ρ_0, grid, 2 / grid.Δz, Zmax())
    @test q[:ρ_0, k] ≈ 2

    q[:ρ_0, 1] = 1.0
    q[:ρ_0, 2] = 2.0
    q[:ρ_0, 3] = 3.0
    q[:a, 1] = 1.0
    q[:a, 2] = 2.0
    q[:a, 3] = 3.0
    ρα_0_cut = q[:ρ_0, Cut(2)] .* q[:a, Cut(2)]
    @test all(ρα_0_cut .== ρα_0_cut)
    ρα_0_dual = q[:ρ_0, Dual(2)] .* q[:a, Dual(2)]
    @test all(ρα_0_dual .== ρα_0_dual)

    q[:a, 3] = 1
    @test !isnan(q)
    @test !isinf(q)
    q[:a, 3] = 1 / 0
    @test isinf(q)
    q[:a, 3] = NaN
    @test isnan(q)

    assign!(q, grid, 2.0)
    vals = [
        q[ϕ, k, i] for ϕ in var_names(q) for i in over_sub_domains(q, ϕ)
        for k in over_elems(grid)
    ]
    @test all(
        [
            q[ϕ, k, i] for ϕ in var_names(q) for i in over_sub_domains(q, ϕ)
            for k in over_elems(grid)
        ] .≈ 2.0,
    )

    assign!(q, grid, 0.0)

    assign!(q, (:w), grid, 2.0)
    @test all(
        [
            q[:w, k, i] for i in over_sub_domains(q, :w)
            for k in over_elems(grid)
        ] .≈ 2.0,
    )

    assign!(q, grid, 0.0)

    assign_real!(q, :w, grid, [2.0 for x in over_elems_real(grid)], gm)
    @test all(
        [
            q[:w, k, gm] for i in over_sub_domains(q, :w)
            for k in over_elems_real(grid)
        ] .≈ 2.0,
    )
    @test all(
        [
            q[:w, k, gm] for i in over_sub_domains(q, :w)
            for k in over_elems_ghost(grid)
        ] .≈ 0.0,
    )

    assign!(q, grid, 2.0)
    assign_ghost!(q, :w, grid, 0.0, gm)
    @test all(
        [
            q[:w, k, gm] for i in over_sub_domains(q, :w)
            for k in over_elems_real(grid)
        ] .≈ 2.0,
    )
    @test all(
        [
            q[:w, k, gm] for i in over_sub_domains(q, :w)
            for k in over_elems_ghost(grid)
        ] .≈ 0.0,
    )

    assign!(q_compare, var_names(q_compare), grid, FT(2.0))
    assign!(q, var_names(q_compare), grid, FT(1.0))

    assign!(q, grid, 0.0)
    assign_real!(q, :w, grid, [2.0 for x in over_elems_real(grid)], gm)
    extrap_0th_order!(q, :w, grid, gm)
    @test all([q[:w, k, gm] for k in over_elems(grid)] .≈ 2.0)

    for k in over_elems(grid)
        q[:a, k] = 2
    end

    assign_ghost!(q, :a, grid, 0.0)
    @test all(q[:a, k] ≈ 0.0 for k in over_elems_ghost(grid))

    extrap!(q, :a, grid)
    k = 1 + grid.n_ghost
    @test q[:a, k - 1] ≈ 2 * q[:a, k] - q[:a, k + 1]
    k = grid.n_elem - grid.n_ghost
    @test q[:a, k + 1] ≈ 2 * q[:a, k] - q[:a, k - 1]

    assign!(q, grid, 0.0)
    ki_1 = first_interior(grid, Zmin())
    ki_2 = first_interior(grid, Zmax())
    kg_1 = first_ghost(grid, Zmin())
    kg_2 = first_ghost(grid, Zmax())

    q[:w, ki_1, gm] = 1.0
    q[:w, ki_2, gm] = 1.0
    extrap!(q, :w, grid, gm)
    @test q[:w, kg_1, gm] ≈ 2.0
    @test q[:w, kg_2, gm] ≈ 2.0

    # Extrapolate field
    q[:w, ki_1, gm] = 2.0
    q[:w, kg_2, gm] = 2.0
    extrap!(q, :w, grid, Zmin(), gm)
    q[:w, ki_1, gm] ≈ 4.0
    q[:w, kg_2, gm] ≈ 2.0

    q[:w, kg_1, gm] = 2.0
    q[:w, ki_2, gm] = 2.0
    extrap!(q, :w, grid, Zmax(), gm)
    q[:w, kg_1, gm] ≈ 2.0
    q[:w, kg_2, gm] ≈ 4.0

    # Extrapolate field with specified boundary value
    q[:w, ki_1, gm] = 0.0
    q[:w, kg_2, gm] = 0.0
    extrap!(q, :w, grid, 1.0, Zmin(), gm)
    q[:w, ki_1, gm] ≈ 2.0
    q[:w, kg_2, gm] ≈ 0.0

    q[:w, kg_1, gm] = 0.0
    q[:w, ki_2, gm] = 0.0
    extrap!(q, :w, grid, 1.0, Zmax(), gm)
    q[:w, kg_1, gm] ≈ 0.0
    q[:w, kg_2, gm] ≈ 2.0

    @test var_string(q, :a, gm) == "a_gm"
    @test var_string(q, :a, en) == "a_en"
    @test var_string(q, :a, ud[1]) == "a_ud_1"
    @test var_string(q, :a, ud[2]) == "a_ud_2"
    @test var_suffix(q, :a, gm) == "_gm"
    @test var_suffix(q, :a, en) == "_en"
    @test var_suffix(q, :a, ud[1]) == "_ud_1"
    @test var_suffix(q, :a, ud[2]) == "_ud_2"
    @test_throws BoundsError var_suffix(q, :a, 1000)
end

output_root = joinpath("..", "output", "tests", "StateVecFuncs")

n_subdomains = 3 # number of sub-domains
n_elems_real = 10 # number of elements

domain_set = DomainSet(gm = 1, en = 1, ud = 1)
grid = UniformGrid(0.0, 1.0, n_elems_real)
vars = (
    (:ρ_0, DomainSubSet(gm = true)),
    (:a, DomainSubSet(gm = true, en = true, ud = true)),
    (:w, DomainSubSet(gm = true, en = true, ud = true)),
    (:ϕ, DomainSubSet(gm = true, en = true, ud = true)),
    (:ψ, DomainSubSet(gm = true, en = true, ud = true)),
)
q = StateVec(vars, grid, domain_set)
vars = (
    (:cv_ϕ_ψ, DomainSubSet(gm = true, en = true, ud = true)),
    (:TCV_ϕ_ψ, DomainSubSet(gm = true)),
)
aux = StateVec(vars, grid, domain_set)

gm, en, ud, sd, al = allcombinations(q)

@testset "IO" begin
    for k in over_elems(grid)
        q[:ρ_0, k] = 1.0 * k
        for i in al
            q[:w, k, i] = 2.0 * k
            q[:a, k, i] = 3.0 * k
        end
    end

    if !Sys.iswindows()
        nc = init_data(NetCDFWriter("q"), grid, q)
        @test isfile("q.nc")

        append_data(nc, q)
        ds = Dataset("q.nc", "r")
        @test all(ds["w_gm"][:] .≈ [q[:w, k, gm] for k in over_elems(grid)])
        @test all(ds["a_gm"][:] .≈ [q[:a, k, gm] for k in over_elems(grid)])
        rm("q.nc")

        export_state(NetCDFWriter("q"), grid, q)
        ds = Dataset("q.nc", "r")
        @test all(ds["w_gm"][:] .≈ [q[:w, k, gm] for k in over_elems(grid)])
        @test all(ds["a_gm"][:] .≈ [q[:a, k, gm] for k in over_elems(grid)])
        rm("q.nc")
    end
end

@testset "Domain average, distribute, diagnose environment, total covariance" begin
    q[:a, 2, en] = 0.25
    q[:a, 2, ud[1]] = 0.75
    q[:w, 2, en] = 2
    q[:w, 2, ud[1]] = 2
    grid_mean!(q, q, :a, :w, grid)
    @test q[:w, 2, gm] ≈ 2

    q[:w, 2, gm] = 2
    distribute!(q, grid, :w)
    @test q[:w, 2, en] ≈ 2
    for i in ud
        @test q[:w, 2, i] ≈ 2
    end

    q[:a, 2, ud[1]] = 0.1
    q[:a, 2, en] = 0.9
    q[:a, 2, gm] = 1.0

    q[:w, 2, gm] = 2
    for i in ud
        q[:w, 2, i] = 2
    end

    diagnose_environment!(q, grid, :a, (:w))
    @test q[:w, 2, en] ≈ 2

    q[:a, 2, 1] = 0.1
    q[:a, 2, 2] = 0.2
    q[:ϕ, 2, 1] = 1
    q[:ϕ, 2, 2] = 2
    q[:ψ, 2, 1] = 2
    q[:ψ, 2, 2] = 3
    aux[:cv_ϕ_ψ, 2, 1] = 1.0
    aux[:cv_ϕ_ψ, 2, 2] = 1.0
    decompose_ϕ_ψ(tcv) = tcv == :TCV_ϕ_ψ ? (:ϕ, :ψ) : error("Bad init")
    total_covariance!(aux, q, aux, :TCV_ϕ_ψ, :cv_ϕ_ψ, :a, grid, decompose_ϕ_ψ)
    @test aux[:TCV_ϕ_ψ, 2] ≈ 0.32
end


@testset "Boundary conditions" begin

    n_subdomains = 3 # number of sub-domains
    n_elems_real = 10 # number of elements

    domain_set = DomainSet(gm = 1, en = 1, ud = 1)
    grid = UniformGrid(0.0, 1.0, n_elems_real)
    FT = eltype(grid)
    vars = (
        (:a, DomainSubSet(gm = true, en = true, ud = true)),
        (:w, DomainSubSet(gm = true, en = true, ud = true)),
    )
    q = StateVec(vars, grid, domain_set)
    vars = (
        (:ρ_0, DomainSubSet(gm = true, en = true, ud = true)),
        (:K, DomainSubSet(gm = true, en = true, ud = true)),
    )
    aux = StateVec(vars, grid, domain_set)

    gm, en, ud, sd, al = allcombinations(q)
    assign!(aux, :K, grid, FT(3))
    assign!(aux, :ρ_0, grid, FT(2))
    assign!(q, :a, grid, FT(0.1))
    @test bc_source(
        q,
        grid,
        aux,
        :w,
        :ρ_0,
        :a,
        :K,
        Zmin(),
        FT(2),
        Dirichlet(),
        DiffusionAbsorbed(),
    ) ≈ FT(240)
    @test bc_source(
        q,
        grid,
        aux,
        :w,
        :ρ_0,
        :a,
        :K,
        Zmin(),
        FT(2),
        Neumann(),
        DiffusionAbsorbed(),
    ) ≈ FT(-4)

    # TODO: Validate bc_source for `AdvectionAbsorbed`
    # @test bc_source(q, grid, aux, :w, :ρ_0, :a, :K,
    #                 Zmin(), FT(2), Dirichlet(), AdvectionAbsorbed()) ≈
    # @test bc_source(q, grid, aux,
    #                 :w, :ρ_0, :a, :K,
    #                 Zmin(), FT(2), Neumann(), AdvectionAbsorbed()) ≈

end
