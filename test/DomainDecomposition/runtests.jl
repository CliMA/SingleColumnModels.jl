using Test
using SingleColumnModels.DomainDecompositions
DD = DomainDecompositions

@testset "DomainSubSet, unit tests" begin

    @test DD.n_subdomains(DD.GridMean{1}) == 1
    @test DD.n_subdomains(DD.Environment{1}) == 1
    @test DD.n_subdomains(DD.Updraft{5}) == 5

    @test DD.gridmean(DomainSubSet(gm = true)) == true
    @test DD.environment(DomainSubSet(gm = true)) == false
    @test DD.updraft(DomainSubSet(gm = true)) == false

    @test DD.gridmean(DomainSubSet(en = true)) == false
    @test DD.environment(DomainSubSet(en = true)) == true
    @test DD.updraft(DomainSubSet(en = true)) == false

    @test DD.gridmean(DomainSubSet(ud = true)) == false
    @test DD.environment(DomainSubSet(ud = true)) == false
    @test DD.updraft(DomainSubSet(ud = true)) == true

    domain_set = DomainSet(gm = 1, en = 1, ud = 4)
    @test DD.n_subdomains(domain_set, DomainSubSet(gm = true)) == (1, 0, 0)
    @test DD.n_subdomains(domain_set, DomainSubSet(en = true)) == (0, 1, 0)
    @test DD.n_subdomains(domain_set, DomainSubSet(ud = true)) == (0, 0, 4)

    @test DD.gridmean(domain_set, DomainSubSet(gm = true)) == 1
    @test DD.environment(domain_set, DomainSubSet(en = true)) == 1
    @test DD.updraft(domain_set, DomainSubSet(ud = true)) == 4

    @test DD.gridmean(domain_set, DomainSubSet(en = true)) == 0
    @test DD.environment(domain_set, DomainSubSet(gm = true)) == 0
    @test DD.updraft(domain_set, DomainSubSet(gm = true)) == 0

    @test DD.sum(domain_set, DomainSubSet(gm = true)) == 1
    @test DD.sum(domain_set, DomainSubSet(en = true)) == 1
    @test DD.sum(domain_set, DomainSubSet(ud = true)) == 4
end

@testset "DomainSubSet, single domain" begin
    domain_set, domain_subset = DomainSet(gm = 1), DomainSubSet(gm = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 1
    @test en == 0
    @test ud == (0,)
    @test idx_ss == idx

    domain_set, domain_subset = DomainSet(en = 1), DomainSubSet(en = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test en == 1
    @test gm == 0
    @test ud == (0,)
    @test idx_ss == idx

    domain_set, domain_subset = DomainSet(ud = 1), DomainSubSet(ud = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test en == 0
    @test gm == 0
    @test ud == (1,)
    @test idx_ss == idx

    domain_set, domain_subset = DomainSet(ud = 2), DomainSubSet(ud = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test en == 0
    @test gm == 0
    @test ud == (1, 2)
    @test idx_ss == idx
end

@testset "DomainSubSet, two domains" begin
    domain_set, domain_subset =
        DomainSet(gm = 1, en = 1), DomainSubSet(gm = true, en = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 2
    @test en == 1
    @test ud == (0,)
    @test idx_ss == idx

    domain_set, domain_subset =
        DomainSet(gm = 1, ud = 3), DomainSubSet(gm = true, ud = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 4
    @test en == 0
    @test ud == (1, 2, 3)
    @test idx_ss == idx

    domain_set, domain_subset =
        DomainSet(en = 1, ud = 3), DomainSubSet(en = true, ud = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 0
    @test en == 4
    @test ud == (1, 2, 3)
    @test idx_ss == idx

    domain_set, domain_subset =
        DomainSet(gm = 1, ud = 1), DomainSubSet(gm = true, ud = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 2
    @test en == 0
    @test ud == (1,)
    @test idx_ss == idx

    domain_set, domain_subset =
        DomainSet(en = 1, ud = 1), DomainSubSet(en = true, ud = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 0
    @test en == 2
    @test ud == (1,)
    @test idx_ss == idx
end

@testset "DomainIdx, all domains" begin
    domain_set, domain_subset = DomainSet(gm = 1, en = 1, ud = 4),
    DomainSubSet(gm = true, en = true, ud = true)
    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 6
    @test en == 5
    @test ud == (1, 2, 3, 4)
    @test idx_ss == idx
end

@testset "DomainIdx, utilizing DomainSubSet" begin
    domain_set = DomainSet(gm = 1, en = 1, ud = 4)
    idx = DomainIdx(domain_set)

    gm, en, ud, sd, al = allcombinations(idx)
    @test gm == 6
    @test en == 5
    @test ud == (1, 2, 3, 4)

    domain_subset = DomainSubSet(gm = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)

    gm, en, ud, sd, al = allcombinations(idx_ss)
    @test gm == 1
    @test en == 0
    @test ud == (0,)

    domain_subset = DomainSubSet(en = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)

    gm, en, ud, sd, al = allcombinations(idx_ss)
    @test gm == 0
    @test en == 1
    @test ud == (0,)

    domain_subset = DomainSubSet(ud = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)

    gm, en, ud, sd, al = allcombinations(idx_ss)
    @test gm == 0
    @test en == 0
    @test ud == (1, 2, 3, 4)
end

@testset "Test DomainSubSet indexing, multiple domains" begin
    domain_set, domain_subset =
        DomainSet(gm = 1, en = 1, ud = 4), DomainSubSet(gm = true)

    idx = DomainIdx(domain_set)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    gm, en, ud, sd, al = allcombinations(idx)

    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    @test DD.updraft(idx_ss) == (0,)
    @test DD.environment(idx_ss) == 0
    @test DD.gridmean(idx_ss) == 1

    domain_subset = DomainSubSet(en = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    @test DD.updraft(idx_ss) == (0,)
    @test DD.environment(idx_ss) == 1
    @test DD.gridmean(idx_ss) == 0

    domain_subset = DomainSubSet(ud = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    @test DD.updraft(idx_ss) == (1, 2, 3, 4)
    @test DD.environment(idx_ss) == 0
    @test DD.gridmean(idx_ss) == 0

    domain_subset = DomainSubSet(gm = true, en = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    @test DD.updraft(idx_ss) == (0,)
    @test DD.environment(idx_ss) == 1
    @test DD.gridmean(idx_ss) == 2

    domain_subset = DomainSubSet(gm = true, ud = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    @test DD.updraft(idx_ss) == (1, 2, 3, 4)
    @test DD.environment(idx_ss) == 0
    @test DD.gridmean(idx_ss) == 5

    domain_subset = DomainSubSet(en = true, ud = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    @test DD.updraft(idx_ss) == (1, 2, 3, 4)
    @test DD.environment(idx_ss) == 5
    @test DD.gridmean(idx_ss) == 0

    domain_subset = DomainSubSet(gm = true, en = true, ud = true)
    idx_ss = DomainIdx(domain_set, domain_subset)
    a_map = DD.get_sv_a_map(idx, idx_ss)
    @test DD.updraft(idx_ss) == (1, 2, 3, 4)
    @test DD.environment(idx_ss) == 5
    @test DD.gridmean(idx_ss) == 6
end

@testset "Test global indexing, many updrafts" begin

    vars = (
        (:ρ_0, DomainSubSet(gm = true)),
        (:a, DomainSubSet(gm = true, en = true, ud = true)),
        (:tke, DomainSubSet(en = true, ud = true)),
        (:K_h, DomainSubSet(ud = true)),
    )

    D = Dict([x => y for (x, y) in vars])

    domain_set = DomainSet(gm = 1, en = 1, ud = 4)
    idx = DomainIdx(domain_set)
    @test DD.n_subdomains(domain_set) == (1, 1, 4)
    @test DD.has_gridmean(idx) == true
    @test DD.has_environment(idx) == true
    @test DD.has_updraft(idx) == true

    sd = subdomains(idx)
    al = alldomains(idx)
    gm, en, ud = eachdomain(idx)
    gm, en, ud, sd, al = allcombinations(idx)

    @test sum(domain_set) == 6

    sd_unmapped = DD.get_sd_unmapped(vars, idx, domain_set)
    sd_mapped = DD.get_sd_mapped(vars, idx, domain_set)

    @test sd_unmapped[:a] == Int[6, 5, 1, 2, 3, 4]
    @test sd_unmapped[:tke] == Int[5, 1, 2, 3, 4]
    @test sd_unmapped[:K_h] == Int[1, 2, 3, 4]
    @test sd_unmapped[:ρ_0] == Int[6]

    @test sd_mapped[:a] == Int[6, 5, 1, 2, 3, 4]
    @test sd_mapped[:tke] == Int[0, 5, 1, 2, 3, 4]
    @test sd_mapped[:K_h] == Int[0, 0, 1, 2, 3, 4]
    @test sd_mapped[:ρ_0] == Int[1, 0, 0, 0, 0]

    vm, dss_per_var, var_names = DD.get_var_mapper(vars, domain_set)

    idx_ss_per_var = Dict([name => DomainIdx(domain_set, dss_per_var[name]) for
    name in var_names]...)
    a_map = Dict([name => DD.get_sv_a_map(idx, idx_ss_per_var[name]) for
    name in var_names]...)
    gm, en, ud, sd, al = allcombinations(idx)

    @test DD.var_string(vm, idx, idx_ss_per_var[:a], :a, gm) == "a_gm"
    @test DD.var_string(vm, idx, idx_ss_per_var[:a], :a, en) == "a_en"
    @test DD.var_string(vm, idx, idx_ss_per_var[:a], :a, ud[1]) == "a_ud_1"
    @test DD.var_string(vm, idx, idx_ss_per_var[:a], :a, ud[2]) == "a_ud_2"
    @test DD.var_string(vm, idx, idx_ss_per_var[:a], :a, ud[3]) == "a_ud_3"
    @test DD.var_string(vm, idx, idx_ss_per_var[:a], :a, ud[4]) == "a_ud_4"
    @test DD.var_suffix(vm, idx, idx_ss_per_var[:a], :a, gm) == "_gm"
    @test DD.var_suffix(vm, idx, idx_ss_per_var[:a], :a, en) == "_en"
    @test DD.var_suffix(vm, idx, idx_ss_per_var[:a], :a, ud[1]) == "_ud_1"
    @test DD.var_suffix(vm, idx, idx_ss_per_var[:a], :a, ud[2]) == "_ud_2"
    @test DD.var_suffix(vm, idx, idx_ss_per_var[:a], :a, ud[3]) == "_ud_3"
    @test DD.var_suffix(vm, idx, idx_ss_per_var[:a], :a, ud[4]) == "_ud_4"
    @test_throws BoundsError DD.var_suffix(
        vm,
        idx,
        idx_ss_per_var[:a],
        :a,
        1000,
    )

    @test DD.get_i_var(a_map[:ρ_0], gm) == 1
    @test DD.get_i_var(a_map[:a], ud[1]) == 1
    @test DD.get_i_var(a_map[:a], ud[2]) == 2
    @test DD.get_i_var(a_map[:a], ud[3]) == 3
    @test DD.get_i_var(a_map[:a], ud[4]) == 4
    @test DD.get_i_var(a_map[:a], en) == 5
    @test DD.get_i_var(a_map[:a], gm) == 6

    @test DD.get_i_state_vec(vm, a_map[:ρ_0], :ρ_0, gm) == 1
    @test DD.get_i_state_vec(vm, a_map[:a], :a, ud[1]) == 2
    @test DD.get_i_state_vec(vm, a_map[:a], :a, ud[2]) == 3
    @test DD.get_i_state_vec(vm, a_map[:a], :a, ud[3]) == 4
    @test DD.get_i_state_vec(vm, a_map[:a], :a, ud[4]) == 5
    @test DD.get_i_state_vec(vm, a_map[:a], :a, en) == 6
    @test DD.get_i_state_vec(vm, a_map[:a], :a, gm) == 7

    @test DD.get_i_state_vec(vm, a_map[:tke], :tke, ud[1]) == 8
    @test DD.get_i_state_vec(vm, a_map[:tke], :tke, ud[2]) == 9
    @test DD.get_i_state_vec(vm, a_map[:tke], :tke, ud[3]) == 10
    @test DD.get_i_state_vec(vm, a_map[:tke], :tke, ud[4]) == 11
    @test DD.get_i_state_vec(vm, a_map[:tke], :tke, en) == 12

    @test DD.get_i_state_vec(vm, a_map[:K_h], :K_h, ud[1]) == 13
    @test DD.get_i_state_vec(vm, a_map[:K_h], :K_h, ud[2]) == 14
    @test DD.get_i_state_vec(vm, a_map[:K_h], :K_h, ud[3]) == 15
    @test DD.get_i_state_vec(vm, a_map[:K_h], :K_h, ud[4]) == 16

end
