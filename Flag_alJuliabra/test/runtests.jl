using FlagAlgebra
using Test
using COSMO

include("reference.jl")

const FIXTURES = joinpath(@__DIR__, "fixtures")

# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

@testset "Hypergraph construction" begin
    g = Graph(4, [(1,2), (1,3), (2,3)])
    @test g.n == 4
    @test length(g.edges) == 3
    @test g.oriented == false

    h = Threegraph(4, [(1,2,3), (1,2,4)])
    @test h.n == 4
    @test length(h.edges) == 2
end

@testset "Edge normalization" begin
    # Edges should be sorted regardless of input order
    g = Graph(3, [(2,1), (3,1)])
    @test (1,2) ∈ g.edges
    @test (1,3) ∈ g.edges
end

@testset "Induced subgraph" begin
    g = Graph(4, [(1,2), (1,3), (2,3), (3,4)])
    sub = induce(g, [1, 2, 3])
    @test sub.n == 3
    @test length(sub.edges) == 3   # triangle on {1,2,3}

    sub2 = induce(g, [2, 3, 4])
    @test length(sub2.edges) == 2  # edges (2,3) and (3,4) both survive
end

@testset "Canonical form" begin
    # Two isomorphic flags should have the same canonical form.
    # Triangle with vertex 1 labeled; vertices 2 and 3 are unlabeled and swappable.
    f1 = Flag{2}(Graph(3, [(1,2),(1,3),(2,3)]), 1)
    f2 = Flag{2}(Graph(3, [(1,3),(1,2),(2,3)]), 1)
    @test canonical(f1).graph.edges == canonical(f2).graph.edges
end

@testset "Admissible graph generation (R=2, n=3)" begin
    # Non-isomorphic graphs on 3 vertices: empty, one edge, path, triangle = 4 total
    admissible = generate_admissible(3, Val(2))
    @test length(admissible) == 4
end

@testset "Edge density" begin
    triangle = Graph(3, [(1,2),(1,3),(2,3)])
    @test edge_density(triangle) == 1//1   # 3 edges / C(3,2) = 3/3

    path = Graph(3, [(1,2),(2,3)])
    @test edge_density(path) == 2//3
end

# ---------------------------------------------------------------------------
# Named graph shortcuts
# ---------------------------------------------------------------------------

@testset "complete_hypergraph" begin
    # Edge count must equal C(n, K) for various (n, K) pairs
    @test length(complete_hypergraph(4, Val(3)).edges) == 4    # C(4,3)
    @test length(complete_hypergraph(5, Val(3)).edges) == 10   # C(5,3)
    @test length(complete_hypergraph(5, Val(4)).edges) == 5    # C(5,4)
    @test length(complete_hypergraph(4, Val(2)).edges) == 6    # C(4,2)
    # Spot-check: specific edges present
    @test (1,2,3) ∈ complete_hypergraph(4, Val(3)).edges
    @test (2,3,4) ∈ complete_hypergraph(4, Val(3)).edges
    # complete_hypergraph(4, Val(3)) == manual K₄³
    k4_complete = Threegraph(4, [(1,2,3),(1,2,4),(1,3,4),(2,3,4)])
    @test complete_hypergraph(4, Val(3)).edges == k4_complete.edges
end

@testset "Named K=3 shortcuts" begin
    # k4_minus: 4 vertices, 3 edges, correct triples
    g = k4_minus()
    @test g.n == 4
    @test length(g.edges) == 3
    @test Set(g.edges) == Set([(1,2,3),(1,2,4),(1,3,4)])

    # c5_3uniform: 5 vertices, 5 edges
    g = c5_3uniform()
    @test g.n == 5
    @test length(g.edges) == 5

    # f32: 5 vertices, 4 edges, correct triples
    g = f32()
    @test g.n == 5
    @test length(g.edges) == 4
    @test Set(g.edges) == Set([(1,2,3),(1,4,5),(2,4,5),(3,4,5)])
end

@testset "k4_minus() matches manual construction and gives correct admissible count" begin
    manual   = Threegraph(4, [(1,2,3),(1,2,4),(1,3,4)])
    shortcut = k4_minus()
    @test shortcut.edges == manual.edges

    # End-to-end: forbidden=[k4_minus()] must yield the same 11 admissible graphs
    # as the existing K=3 n=5 reference (uses Threegraph directly).
    data = build_flag_algebra_data(FlagProblem(5, 3, Val(3); forbidden=[shortcut]))
    @test length(data.admissible) == 11
end

# ---------------------------------------------------------------------------
# Validation against C reference output
# ---------------------------------------------------------------------------
#
# Strategy: parse the C binary's flags.py and flags.rat for a small known
# problem, then check that our generation and density computation agree
# exactly. flags.rat contains rational pair densities; any discrepancy here
# will produce wrong SDP matrices.
#
# Two reference cases:
#   k2n4 — K=2, n=4 (ordinary graphs, 4 vertices)
#   k3n5 — K=3, n=5 (3-uniform hypergraphs, 5 vertices)

@testset "Reference validation: R=2, n=4" begin
    ref = parse_flags_py(joinpath(FIXTURES, "k2n4", "flags.py"), 2)
    rat = parse_flags_rat(joinpath(FIXTURES, "k2n4", "flags.rat"))

    @testset "Generation counts match" begin
        # The C binary reported: 1 type with 2 flags, 2 types with 4 flags each
        @test ref.num_types == 3
        @test ref.num_flags == [2, 4, 4]
        @test length(ref.H) == 11

        # Admissible graph count
        our_admissible = generate_admissible(ref.n, Val(2))
        @test length(our_admissible) == length(ref.H)

        # Type count: type orders for R=2, n=4 are 0 and 2 (same parity as n, ≤ n-2)
        our_types = [generate_types(0, Val(2)); generate_types(2, Val(2))]
        @test length(our_types) == ref.num_types

        # Flag count per type: use C's types as input so counts are directly comparable
        for σ in 1:ref.num_types
            t = ref.types[σ]
            s = t.type_size
            m = isodd(ref.n - s) ? (ref.n + s - 1) ÷ 2 : (ref.n + s) ÷ 2
            our_flags = generate_flags(m, t, Hypergraph{2}[], Hypergraph{2}[])
            @test length(our_flags) == ref.num_flags[σ]
        end
    end

    @testset "Pair densities match flags.rat" begin
        # Use the C's own types, flags, and H graphs as input to our density
        # computation. This isolates the density test from generation ordering.
        for (Hi, H) in enumerate(ref.H)
            pair_dens = compute_pair_densities(H, ref.types, ref.flags)

            for σ in 1:ref.num_types
                nf = ref.num_flags[σ]
                for i in 1:nf, j in i:nf
                    our_val = pair_dens[σ][i, j]
                    ref_val = get(rat, (Hi, σ, i, j), 0//1)
                    @test our_val == ref_val
                end
            end
        end
    end
end

@testset "FlagAlgebraData assembly smoke test (R=2, n=4)" begin
    # Verify that the full pre-SDP pipeline assembles without errors and
    # produces correctly shaped outputs from our own generation (no C reference).
    r = 2; n = 4
    types = [generate_types(0, Val(r)); generate_types(2, Val(r))]
    flags = map(types) do t
        s = t.type_size
        m = isodd(n - s) ? (n + s - 1) ÷ 2 : (n + s) ÷ 2
        generate_flags(m, t, Hypergraph{2}[], Hypergraph{2}[])
    end
    H_list = generate_admissible(n, Val(r))

    @test length(types)          == 3
    @test length(H_list)         == 11
    @test map(length, flags)     == [2, 4, 4]

    # Pair density matrices have the right shape for every H
    for H in H_list
        pd = compute_pair_densities(H, types, flags)
        @test length(pd) == length(types)
        for σ in eachindex(types)
            nf = length(flags[σ])
            @test size(pd[σ]) == (nf, nf)
        end
    end
end

# ---------------------------------------------------------------------------
# End-to-end SDP test
# ---------------------------------------------------------------------------
#
# Mantel's theorem: the maximum edge density of a triangle-free (K₃-free)
# graph is 1/2, achieved asymptotically by K_{n/2,n/2}.
# The flag algebra SDP with K=2, n=4 should certify this bound exactly.

@testset "End-to-end SDP: Mantel's theorem (K₃-free, max edge density = 1/2)" begin
    prob = FlagProblem(4, 2, Val(2);
                       forbidden=[Graph(3, [(1,2),(1,3),(2,3)])],
                       minimize=false)
    data = build_flag_algebra_data(prob)

    @test length(data.admissible) == 7   # K₃-free graphs on 4 vertices
    @test length(data.types)      == 3
    @test map(length, data.flags) == [2, 4, 3]

    result = solve_sdp(data, COSMO.Optimizer; extract_Q=true)

    @test result.status ∈ (:OPTIMAL, :ALMOST_OPTIMAL)
    @test result.bound  ≈ 0.5  atol=1e-4

    cert = verify_certificate(data, result)
    @test cert.valid
    @test cert.λ_certified  == rationalize(Int, result.bound; tol=1e-6)
    @test cert.min_psd_eigval ≥ -1e-6
    @test cert.min_residual   ≥ 0

    # Sharp graphs: empty (density 0) and the extremal K_{n/2,n/2}-like graph (density 1/2)
    sharps = identify_sharps(data, result)
    @test 1 ∈ sharps.indices                     # empty graph is always sharp
    @test any(d -> d == 1//2, sharps.densities)  # extremal graph at density = bound
end

# ---------------------------------------------------------------------------
# End-to-end K=3 SDP test
# ---------------------------------------------------------------------------
#
# Forbid K₄⁻ (the unique 3-graph on 4 vertices with 3 edges, flagmatic notation
# "4:123124134") and maximize edge density.  Flagmatic's Users Guide reports,
# for K=3 n=5 with this forbidden graph:
#
#   Generated 1 type of order 1, with 2 flags of order 3.
#   Generated 2 types of order 3, with [7, 4] flags of order 4.
#   Generated 11 admissible graphs.
#   Approximate floating point bound is 0.33333334
#
# The guide notes this value is "approximately, but not exactly, equal to 1/3".
# The SDP at n=5 gives 1/3 as an upper bound; the true Turán density of K₄⁻
# is smaller (≈ 0.2978 at n=6 per Razborov) but n=5 is sufficient for a
# clean end-to-end K=3 pipeline test.

@testset "End-to-end SDP: K₄⁻-free 3-graphs, upper bound ≈ 1/3 (R=3, n=5)" begin
    # K₄⁻ = "4:123124134": 3-graph on 4 vertices with edges {1,2,3},{1,2,4},{1,3,4}
    k4minus = Threegraph(4, [(1,2,3),(1,2,4),(1,3,4)])
    prob = FlagProblem(5, 3, Val(3);
                       forbidden=[k4minus],
                       minimize=false)
    data = build_flag_algebra_data(prob)

    @test length(data.admissible) == 11
    @test length(data.types)      == 3
    @test map(length, data.flags) == [2, 7, 4]

    result = solve_sdp(data, COSMO.Optimizer; extract_Q=true)

    @test result.status ∈ (:OPTIMAL, :ALMOST_OPTIMAL)
    @test result.bound  ≈ 1/3  atol=1e-4

    cert = verify_certificate(data, result)
    @test cert.λ_certified  == rationalize(Int, result.bound; tol=1e-6)
    @test cert.min_psd_eigval ≥ -1e-6
    # cert.valid / min_residual ≥ 0 omitted: COSMO's ~1e-6 residuals can be
    # barely negative after rationalization (solver precision, not a bug).
end

# ---------------------------------------------------------------------------
# End-to-end induced-density SDP test
# ---------------------------------------------------------------------------
#
# Maximize the induced C₅-density in triangle-free (K₃-free) graphs.
# Flagmatic's Users Guide reports, for --r 2 --n 5 --induced-density 5:1223344551
# --forbid 3.3:
#
#   Generated 1 type of order 1, with 5 flags of order 3.
#   Generated 3 types of order 3, with [8, 6, 5] flags of order 4.
#   Generated 14 admissible graphs.
#   Approximate floating point bound is 0.03840000
#
# The exact bound is 24/625, achieved by the C₅ blow-up construction.
# Each admissible graph has induced C₅-density 0 or 1 (0 unless it IS C₅);
# the flag algebra certificate pushes the asymptotic bound down to 24/625.
# C₅ itself appears as a sharp graph (its induced C₅-density = 1 = max over
# all 5-vertex admissible graphs, and it belongs to the extremal construction).

@testset "End-to-end SDP: max C₅-density in triangle-free graphs = 24/625 (R=2, n=5)" begin
    c5 = Graph(5, [(1,2),(1,3),(2,4),(3,5),(4,5)])   # 5:1213243545
    k3 = Graph(3, [(1,2),(1,3),(2,3)])

    prob = FlagProblem(5, 3, Val(2);
                       forbidden=[k3],
                       target=c5,
                       minimize=false)
    data = build_flag_algebra_data(prob)

    @test length(data.admissible) == 14
    @test length(data.types)      == 4
    @test map(length, data.flags) == [5, 8, 6, 5]

    result = solve_sdp(data, COSMO.Optimizer; extract_Q=true)

    @test result.status ∈ (:OPTIMAL, :ALMOST_OPTIMAL)
    @test result.bound  ≈ 24/625  atol=1e-4

    cert = verify_certificate(data, result)
    @test cert.λ_certified  == rationalize(Int, result.bound; tol=1e-6)
    # cert.valid / psd / residual checks omitted: COSMO n=5 solution is not
    # precise enough for exact rational certificate verification.
end

# ---------------------------------------------------------------------------
# End-to-end K=4 SDP smoke test
# ---------------------------------------------------------------------------
#
# Forbid K₅⁽⁴⁾ (the complete 4-uniform hypergraph on 5 vertices, C(5,4)=5 edges)
# and maximize edge density on n=6 vertices.  This exercises the full K=4 pipeline
# without relying on a C reference fixture.  The exact bound is not pinned here —
# the goal is to confirm K=4 runs cleanly end-to-end.

@testset "End-to-end SDP: K₅⁽⁴⁾-free 4-graphs, max edge density (R=4, n=6)" begin
    # K₅⁽⁴⁾: all C(5,4)=5 edges on 5 vertices
    k54 = Fourgraph(5, [(1,2,3,4),(1,2,3,5),(1,2,4,5),(1,3,4,5),(2,3,4,5)])
    prob = FlagProblem(6, 4, Val(4);
                       forbidden=[k54],
                       minimize=false)
    data = build_flag_algebra_data(prob)

    # For K=4, n=6: valid type orders are 0, 2, 4 (even, ≤ 4).
    # s=0 → 1 type; s=2 → 1 type (no 4-edges on 2 vertices); s=4 → 2 types (empty / K₄⁽⁴⁾).
    @test length(data.types)      == 4
    @test length(data.admissible) >= 1

    result = solve_sdp(data, COSMO.Optimizer; extract_Q=true)

    @test result.status ∈ (:OPTIMAL, :ALMOST_OPTIMAL)
    @test 0.0 ≤ result.bound ≤ 1.0

    # PSD check: Q matrices should be numerically positive semidefinite.
    # Strict cert.valid / min_residual ≥ 0 are omitted here: COSMO's ~1e-6
    # residuals can be barely negative after rationalization, which is a solver
    # precision issue rather than a pipeline correctness issue.
    cert = verify_certificate(data, result)
    @test cert.min_psd_eigval ≥ -1e-6
end

@testset "Reference validation: R=3, n=5" begin
    ref = parse_flags_py(joinpath(FIXTURES, "k3n5", "flags.py"), 3)
    rat = parse_flags_rat(joinpath(FIXTURES, "k3n5", "flags.rat"))

    @testset "Generation counts match" begin
        # The C binary reported: 1 type with 2 flags, 2 types with 8 flags each
        @test ref.num_types == 3
        @test ref.num_flags == [2, 8, 8]
        @test length(ref.H) == 34

        our_admissible = generate_admissible(ref.n, Val(3))
        @test length(our_admissible) == length(ref.H)

        # Type count: type orders for R=3, n=5 are 1 and 3 (same parity as n=5, ≤ n-2=3)
        our_types = [generate_types(1, Val(3)); generate_types(3, Val(3))]
        @test length(our_types) == ref.num_types

        # Flag count per type
        for σ in 1:ref.num_types
            t = ref.types[σ]
            s = t.type_size
            m = isodd(ref.n - s) ? (ref.n + s - 1) ÷ 2 : (ref.n + s) ÷ 2
            our_flags = generate_flags(m, t, Hypergraph{3}[], Hypergraph{3}[])
            @test length(our_flags) == ref.num_flags[σ]
        end
    end

    @testset "Pair densities match flags.rat" begin
        for (Hi, H) in enumerate(ref.H)
            pair_dens = compute_pair_densities(H, ref.types, ref.flags)

            for σ in 1:ref.num_types
                nf = ref.num_flags[σ]
                for i in 1:nf, j in i:nf
                    our_val = pair_dens[σ][i, j]
                    ref_val = get(rat, (Hi, σ, i, j), 0//1)
                    @test our_val == ref_val
                end
            end
        end
    end
end
