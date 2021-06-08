heat_eq_dir = joinpath(output_root, "HeatEquation")

@testset "∂_t T = K ΔT + 1, T = 0 ∈ ∂Ω, diffusion-explicit" begin
    domain_set = DomainSet(gm = 1)
    domain_subset = DomainSubSet(gm = true)
    K, maxiter, Δt = 1.0, 1000, 0.005
    grid = UniformGrid(0.0, 1.0, 10)
    q = StateVec(((:T, domain_subset),), grid, domain_set)
    tmp = StateVec(((:ΔT, domain_subset),), grid, domain_set)
    rhs = deepcopy(q)
    for i in 1:maxiter
        for k in over_elems_real(grid)
            tmp[:ΔT, k] = Δ_z(q[:T, Cut(k)], grid)
            rhs[:T, k] = K * tmp[:ΔT, k] + 1
        end
        for k in over_elems(grid)
            q[:T, k] += Δt * rhs[:T, k]
        end
        apply_Dirichlet!(q, :T, grid, 0.0, Zmax())
        apply_Dirichlet!(q, :T, grid, 0.0, Zmin())
    end
    sol_analtyic = 1 / 2 * grid.zc - grid.zc .^ 2 / 2
    sol_error = [abs(q[:T, k] - sol_analtyic[k]) for k in over_elems_real(grid)]
    @test all(sol_error .< grid.Δz)
    mkpath(heat_eq_dir)
    nc = NetCDFWriter(joinpath(heat_eq_dir, "T"))
    export_state(nc, grid, q)
    rm(full_name(nc))
end


@testset "∂_t T = K ΔT, T={1, 0} ∈ {z_{min}, z_{max}}, diffusion-explicit" begin
    domain_set = DomainSet(gm = 1)
    domain_subset = DomainSubSet(gm = true)
    K, maxiter, Δt = 1.0, 1000, 0.002
    grid = UniformGrid(0.0, 1.0, 10)
    q = StateVec(((:T, domain_subset),), grid, domain_set)
    tmp = StateVec(((:ΔT, domain_subset),), grid, domain_set)
    rhs = deepcopy(q)
    for i in 1:maxiter
        for k in over_elems_real(grid)
            tmp[:ΔT, k] = Δ_z(q[:T, Cut(k)], grid)
            rhs[:T, k] = K * tmp[:ΔT, k]
        end
        for k in over_elems(grid)
            q[:T, k] += Δt * rhs[:T, k]
        end
        apply_Dirichlet!(q, :T, grid, 1.0, Zmin())
        apply_Dirichlet!(q, :T, grid, 0.0, Zmax())
    end
    sol_analtyic = 1 .- grid.zc
    sol_error = [abs(q[:T, k] - sol_analtyic[k]) for k in over_elems_real(grid)]
    @test all(sol_error .< grid.Δz)
end


@testset "∂_t T = K ΔT, T = {1, 0} ∈ {z_{min}, z_{max}}, diffusion implicit" begin
    domain_set = DomainSet(gm = 1)
    domain_subset = DomainSubSet(gm = true)
    K, maxiter, Δt = 1.0, 1000, 0.1
    grid = UniformGrid(0.0, 1.0, 10)
    k_star1 = first_interior(grid, Zmin())
    k_star2 = first_interior(grid, Zmax())
    vars = ((:T, domain_subset), (:a, domain_subset))
    q = StateVec(vars, grid, domain_set)
    vars = (
        (:ΔT, domain_subset),
        (:a_tau, domain_subset),
        (:ρ_0, domain_subset),
        (:K, domain_subset),
        (:a, domain_subset),
    )
    tmp = StateVec(vars, grid, domain_set)
    rhs = deepcopy(q)
    for k in over_elems(grid)
        q[:a, k] = 1.0
        tmp[:ρ_0, k] = 1.0
        tmp[:a_tau, k] = 1.0
        tmp[:K, k] = 1.0
    end
    for i in 1:maxiter
        for k in over_elems_real(grid)
            rhs[:T, k] = 0.0
            if k == k_star1
                gv = 2 * 1 - q[:T, k]
                grad_T = ∇_z_flux([gv, 0], grid)
                bc_surf = ∇_z_flux([grad_T, 0], grid)
                rhs[:T, k] += bc_surf
            end
            if k == k_star2
                gv = -q[:T, k]
                grad_T = ∇_z_flux([0, gv], grid)
                bc_surf = ∇_z_flux([0, grad_T], grid)
                rhs[:T, k] += bc_surf
            end
        end
        solve_tdma!(q, rhs, tmp, :T, :ρ_0, :a, :a_tau, :K, grid, Δt)
        extrap!(q, :T, grid)
    end
    sol_analtyic = 1 .- grid.zc
    sol_numerical = [q[:T, k] for k in over_elems_real(grid)]
    sol_error = [abs(q[:T, k] - sol_analtyic[k]) for k in over_elems_real(grid)]
    @test all(sol_error .< grid.Δz)
end


@testset "∂_t T = K ΔT, -K∇T = F ∈ z_{min}, T = 0 ∈ z_{max}, explicit Euler" begin
    domain_set = DomainSet(gm = 1)
    domain_subset = DomainSubSet(gm = true)
    K, maxiter, Δt = 1.0, 10, 0.001, 10
    grid = UniformGrid(0.0, 1.0, 10)
    k_star1 = first_interior(grid, Zmin())
    k_star2 = first_interior(grid, Zmax())
    q = StateVec(((:T1, domain_subset), (:T2, domain_subset)), grid, domain_set)
    tmp = StateVec(
        ((:ΔT1, domain_subset), (:ΔT2, domain_subset)),
        grid,
        domain_set,
    )
    rhs = deepcopy(q)
    for i in 1:maxiter
        q_out = -K * 1
        q_in = -q_out
        for k in over_elems_real(grid)
            tmp[:ΔT1, k] = Δ_z(q[:T1, Cut(k)], grid)
            tmp[:ΔT2, k] = Δ_z(q[:T2, Cut(k)], grid)
            rhs[:T1, k] = K * tmp[:ΔT1, k]
            rhs[:T2, k] = K * tmp[:ΔT2, k]
            if k == k_star1 && i > 1
                rhs[:T2, k] += -∇_z_flux([q_in, 0.0], grid)
            end
        end
        for k in over_elems(grid)
            q[:T1, k] += Δt * rhs[:T1, k]
            q[:T2, k] += Δt * rhs[:T2, k]
        end
        apply_Neumann!(q, :T1, grid, q_out, Zmin())
        apply_Dirichlet!(q, :T1, grid, 0.0, Zmax())
        apply_Neumann!(q, :T2, grid, 0.0, Zmin())
        apply_Dirichlet!(q, :T2, grid, 0.0, Zmax())
        @test all([q[:T1, k] ≈ q[:T2, k] for k in over_elems_real(grid)])
    end
end


@testset "∂_t T = ∇ • (K(z) ∇T) + 1, T = 0 ∈ ∂Ω, explicit Euler" begin
    domain_set = DomainSet(gm = 1)
    domain_subset = DomainSubSet(gm = true)
    K, maxiter, Δt = 1.0, 1000, 0.001
    grid = UniformGrid(0.0, 1.0, 10)
    unknowns = ((:T, domain_subset),)
    vars = ((:ΔT, domain_subset), (:K_thermal, domain_subset))
    q = StateVec(unknowns, grid, domain_set)
    tmp = StateVec(vars, grid, domain_set)
    rhs = deepcopy(q)
    cond_thermal(z) = z > 0.5 ? 1 : 0.1
    for k in over_elems(grid)
        tmp[:K_thermal, k] = cond_thermal(grid.zc[k])
    end
    for i in 1:maxiter
        for k in over_elems_real(grid)
            tmp[:ΔT, k] = Δ_z(q[:T, Cut(k)], grid, tmp[:K_thermal, Cut(k)])
            rhs[:T, k] = K * tmp[:ΔT, k] + 1
        end
        for k in over_elems(grid)
            q[:T, k] += Δt * rhs[:T, k]
        end
        apply_Dirichlet!(q, :T, grid, 0.0, Zmax())
        apply_Dirichlet!(q, :T, grid, 0.0, Zmin())
    end
end


@testset "∂_t T = ∇ • (K(z) ∇T) + 1, T = {0, 0} ∈ {z_{min}, z_{max}}, comparison" begin
    domain_set = DomainSet(gm = 1)
    domain_subset = DomainSubSet(gm = true)
    VS, TS, F, K = 1.0, 1.0, 1.0, 1.0
    maxiter, Δt = 5000, 0.005
    grid = UniformGrid(0.0, 1.0, 10)
    k_star1 = first_interior(grid, Zmin())
    k_star2 = first_interior(grid, Zmax())
    vars = (
        (:T_explicit_surf, domain_subset),
        (:T_explicit_vol, domain_subset),
        (:T_implicit, domain_subset),
        (:a, domain_subset),
    )
    q = StateVec(vars, grid, domain_set)
    vars = (
        (:ΔT, domain_subset),
        (:a_tau, domain_subset),
        (:ρ_0, domain_subset),
        (:K_thermal, domain_subset),
        (:a, domain_subset),
    )
    tmp = StateVec(vars, grid, domain_set)
    rhs = deepcopy(q)
    for k in over_elems(grid)
        q[:a, k] = 1.0
        tmp[:ρ_0, k] = 1.0
        tmp[:a_tau, k] = 1.0
    end
    cond_thermal(z) = z > 0.5 ? 1 : 0.1
    for k in over_elems(grid)
        tmp[:K_thermal, k] = cond_thermal(grid.zc[k])
    end
    for i in 1:maxiter
        for k in over_elems_real(grid)
            tmp[:ΔT, k] =
                Δ_z(q[:T_explicit_surf, Cut(k)], grid, tmp[:K_thermal, Cut(k)])
            rhs[:T_explicit_surf, k] = tmp[:ΔT, k] + VS
            tmp[:ΔT, k] =
                Δ_z(q[:T_explicit_vol, Cut(k)], grid, tmp[:K_thermal, Cut(k)])
            rhs[:T_explicit_vol, k] = tmp[:ΔT, k] + VS
            rhs[:T_implicit, k] = VS
            k == k_star1 && (
                rhs[:T_implicit, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_implicit,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmin(),
                    TS,
                    Dirichlet(),
                    DiffusionAbsorbed(),
                )
            )
            k == k_star2 && (
                rhs[:T_implicit, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_implicit,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmax(),
                    TS,
                    Dirichlet(),
                    DiffusionAbsorbed(),
                )
            )
            k == k_star1 && (
                rhs[:T_explicit_vol, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_explicit_vol,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmin(),
                    TS,
                    Dirichlet(),
                    DiffusionAbsorbed(),
                )
            )
            k == k_star2 && (
                rhs[:T_explicit_vol, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_explicit_vol,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmax(),
                    TS,
                    Dirichlet(),
                    DiffusionAbsorbed(),
                )
            )
        end

        for k in over_elems_real(grid)
            q[:T_explicit_surf, k] += Δt * rhs[:T_explicit_surf, k]
            q[:T_explicit_vol, k] += Δt * rhs[:T_explicit_vol, k]
        end
        solve_tdma!(
            q,
            rhs,
            tmp,
            :T_implicit,
            :ρ_0,
            :a,
            :a_tau,
            :K_thermal,
            grid,
            Δt,
        )

        apply_Dirichlet!(q, :T_explicit_surf, grid, TS, Zmax())
        apply_Dirichlet!(q, :T_explicit_surf, grid, TS, Zmin())

        assign_ghost!(q, :T_explicit_vol, grid, 0.0, Zmin())
        assign_ghost!(q, :T_explicit_vol, grid, 0.0, Zmax())

        assign_ghost!(q, :T_implicit, grid, 0.0, Zmin())
        assign_ghost!(q, :T_implicit, grid, 0.0, Zmax())
    end
    extrap!(q, :T_explicit_vol, grid, TS, Zmin())
    extrap!(q, :T_explicit_vol, grid, TS, Zmax())
    extrap!(q, :T_implicit, grid, TS, Zmin())
    extrap!(q, :T_implicit, grid, TS, Zmax())

    @test q[:T_explicit_vol, Dual(k_star1)][1] ≈ TS
    @test q[:T_explicit_vol, Dual(k_star2)][2] ≈ TS
    @test all([
        q[:T_explicit_vol, k] ≈ q[:T_explicit_surf, k]
        for k in over_elems_real(grid)
    ])
    @test all([
        q[:T_implicit, k] ≈ q[:T_explicit_surf, k]
        for k in over_elems_real(grid)
    ])

end


@testset "∂_t T = ∇ • (K(z) ∇T) + VS, -K(z) ∇T = 1 ∈ z_{min}, T = 0 ∈ z_{max}, comparison" begin
    domain_set = DomainSet(gm = 1)
    domain_subset = DomainSubSet(gm = true)
    VS, TS, F = 1.0, 0.0, 2.0
    IC = 0.0
    maxiter, Δt = 10000, 0.0001
    grid = UniformGrid(0.0, 1.0, 20)
    k_star1 = first_interior(grid, Zmin())
    k_star2 = first_interior(grid, Zmax())
    vars = (
        (:T_explicit_surf, domain_subset),
        (:T_explicit_vol, domain_subset),
        (:T_implicit, domain_subset),
        (:a, domain_subset),
    )
    q = StateVec(vars, grid, domain_set)
    vars = (
        (:ΔT, domain_subset),
        (:a_tau, domain_subset),
        (:ρ_0, domain_subset),
        (:K_thermal, domain_subset),
    )
    tmp = StateVec(vars, grid, domain_set)
    rhs = deepcopy(q)
    for k in over_elems(grid)
        q[:a, k] = 1.0
        tmp[:ρ_0, k] = 1.0
        tmp[:a_tau, k] = 1.0
    end
    cond_thermal(z) = 5 + 3 * cos(6 * π * z)
    for k in over_elems(grid)
        tmp[:K_thermal, k] = cond_thermal(grid.zc[k])
        q[:T_explicit_surf, k] = IC
        q[:T_explicit_vol, k] = IC
        q[:T_implicit, k] = IC
    end
    K_surf = tmp[:K_thermal, Dual(k_star1)][1]
    apply_Neumann!(q, :T_explicit_surf, grid, F / K_surf, Zmin())
    apply_Dirichlet!(q, :T_explicit_surf, grid, TS, Zmax())
    for i in 1:maxiter
        assign_ghost!(q, :T_explicit_vol, grid, 0.0, Zmin())
        assign_ghost!(q, :T_explicit_vol, grid, 0.0, Zmax())
        for k in over_elems_real(grid)
            tmp[:ΔT, k] = Δ_z_dual(
                q[:T_explicit_surf, Cut(k)],
                grid,
                tmp[:K_thermal, Dual(k)],
            )
            rhs[:T_explicit_surf, k] = tmp[:ΔT, k] + VS

            tmp[:ΔT, k] = Δ_z_dual(
                q[:T_explicit_vol, Cut(k)],
                grid,
                tmp[:K_thermal, Dual(k)],
            )
            rhs[:T_explicit_vol, k] = tmp[:ΔT, k] + VS
            k == k_star1 && (
                rhs[:T_explicit_vol, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_explicit_vol,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmin(),
                    F,
                    Neumann(),
                    DiffusionAbsorbed(),
                )
            )
            k == k_star2 && (
                rhs[:T_explicit_vol, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_explicit_vol,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmax(),
                    TS,
                    Dirichlet(),
                    DiffusionAbsorbed(),
                )
            )

            rhs[:T_implicit, k] = VS

            k == k_star1 && (
                rhs[:T_implicit, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_implicit,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmin(),
                    F,
                    Neumann(),
                    DiffusionAbsorbed(),
                )
            )
            k == k_star2 && (
                rhs[:T_implicit, k] += bc_source(
                    q,
                    grid,
                    tmp,
                    :T_implicit,
                    :ρ_0,
                    :a,
                    :K_thermal,
                    Zmax(),
                    TS,
                    Dirichlet(),
                    DiffusionAbsorbed(),
                )
            )
        end

        for k in over_elems_real(grid)
            q[:T_explicit_surf, k] += Δt * rhs[:T_explicit_surf, k]
            q[:T_explicit_vol, k] += Δt * rhs[:T_explicit_vol, k]
        end
        solve_tdma!(
            q,
            rhs,
            tmp,
            :T_implicit,
            :ρ_0,
            :a,
            :a_tau,
            :K_thermal,
            grid,
            Δt,
        )

        apply_Neumann!(q, :T_explicit_surf, grid, F / K_surf, Zmin())
        apply_Dirichlet!(q, :T_explicit_surf, grid, TS, Zmax())

        apply_Neumann!(q, :T_explicit_vol, grid, F / K_surf, Zmin())
        apply_Dirichlet!(q, :T_explicit_vol, grid, TS, Zmax())

        apply_Neumann!(q, :T_implicit, grid, F / K_surf, Zmin())
        apply_Dirichlet!(q, :T_implicit, grid, TS, Zmax())
    end

    apply_Neumann!(q, :T_explicit_vol, grid, F / K_surf, Zmin())
    extrap!(q, :T_explicit_vol, grid, TS, Zmax())

    apply_Neumann!(q, :T_implicit, grid, F / K_surf, Zmin())
    extrap!(q, :T_implicit, grid, TS, Zmax())

    @test q[:T_explicit_vol, Dual(k_star2)][2] ≈ TS
    @test q[:T_implicit, Dual(k_star2)][2] ≈ TS

    @test all([
        abs(q[:T_explicit_vol, k] - q[:T_explicit_surf, k]) < 10 * Δt
        for k in over_elems(grid)
    ])
    @test all([
        abs(q[:T_explicit_vol, k] - q[:T_explicit_surf, k]) < 10 * Δt
        for k in over_elems(grid)
    ])
    @test all([
        abs(q[:T_implicit, k] - q[:T_explicit_surf, k]) < 10 * Δt
        for k in over_elems(grid)
    ])
end
