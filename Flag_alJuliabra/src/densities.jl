# Density computations: edge density, induced subgraph density, pair densities.
#
# NOTE: pair density normalization
# C's tcnum = P(n,s) × C(n-s,half) × C(n-s-half,half) / 2 — a pure
# combinatorial constant counting ALL (type_perm, half1, half2) combinations,
# NOT restricted to type_perms that induce the correct type graph.
# This implementation uses ordered iteration (both orderings of each half-pair),
# so total = 2 × tcnum. All counts are also doubled relative to C's.
# Diagonal:     counts[i,i] // total        (2x / 2x = same as C's factor-2 diagonal)
# Off-diagonal: counts[i,j] // (2 * total)  (2x / 4x = same as C's c/(2*tcnum))
# See NOTES.md.

using Combinatorics

# Edge density of a hypergraph: m / C(n, R)
function edge_density(g::Hypergraph{R}) :: Rational{Int} where R
    length(g.edges) // binomial(g.n, R)
end

# Density of h as an induced subgraph of g:
# fraction of R-subsets of V(g) whose induced subgraph is isomorphic to h.
function induced_density(g::Hypergraph{R}, h::Hypergraph{R}) :: Rational{Int} where R
    h_canon = canonical(Flag{R}(h, 0)).graph.edges
    num_found = count(combinations(1:g.n, h.n)) do verts
        canonical(Flag{R}(induce(g, collect(verts)), 0)).graph.edges == h_canon
    end
    num_found // binomial(g.n, h.n)
end

# Compute pair density matrices for one admissible graph H.
#
# Returns a vector of upper-triangular matrices (one per type σ).
# pair_dens[σ][i,j] (i ≤ j) = probability that a random injection of the
# type vertices, extended by two independent random flag halves, yields
# flags i and j over type σ.
#
# m = ⌊(n+s)/2⌋ is the flag vertex count, chosen so that two flags
# together cover (or nearly cover) all of H's vertices.
# When n-s is odd, one vertex of H is unused in each pairing.

function compute_pair_densities(H::Hypergraph{R},
                                 types::Vector{Flag{R}},
                                 flags::Vector{Vector{Flag{R}}}) :: Vector{Matrix{Rational{Int}}} where R
    n         = H.n
    num_types = length(types)

    # Precompute canonical forms once to avoid repeating in inner loops.
    type_canons = [canonical(t).graph.edges for t in types]
    flag_canons = [[canonical(f).graph.edges for f in flags[σ]] for σ in 1:num_types]

    map(1:num_types) do σ
        s    = types[σ].type_size
        m    = isodd(n - s) ? (n + s - 1) ÷ 2 : (n + s) ÷ 2
        half = m - s
        nf   = length(flags[σ])

        # Each type has its own counts and total. Types with different s values
        # produce different iteration counts; mixing them in a shared total was
        # the original bug — it contaminated every type's denominator.
        # Types sharing the same s get the same total regardless, since total
        # counts all iterations unconditionally (matching the C's tcnum).
        counts = zeros(Int, nf, nf)
        total  = 0

        for type_perm in permutations(1:n, s)
            # Pre-check type match once per type_perm (hoisted out of inner loops)
            type_sub     = induce(H, collect(type_perm))
            type_matches = type_sub.edges == type_canons[σ]

            rest = setdiff(1:n, type_perm)

            for half1 in combinations(rest, half)
                rest2 = setdiff(rest, half1)
                length(rest2) < half && continue

                for half2 in combinations(rest2, half)
                    # Count ALL (type_perm, half1, half2) regardless of type match.
                    # This mirrors C's tcnum, which is a pure combinatorial constant
                    # P(n,s) × C(n-s,half) × C(n-s-half,half) / 2 (unordered), not
                    # restricted to type_perms that induce the correct type graph.
                    total += 1

                    type_matches || continue

                    f1 = find_flag(H, type_perm, half1, flag_canons[σ], s)
                    f2 = find_flag(H, type_perm, half2, flag_canons[σ], s)
                    (f1 === nothing || f2 === nothing) && continue

                    i, j = minmax(f1, f2)
                    counts[i, j] += 1
                end
            end
        end

        # Normalize counts into rational pair densities.
        #
        # With ordered iteration (both (half1,half2) and (half2,half1) visited):
        #   - total  = 2 × C's tcnum  (every unordered pair counted twice)
        #   - diagonal counts  = 2 × C's counts  (both orderings give same flag pair)
        #   - off-diagonal counts = 2 × C's counts  (same doubling)
        #
        # C normalizes off-diagonal by 2*tcnum and applies factor=2 only for diagonal.
        # In our ordered scheme:
        #   diagonal:     counts[i,i] / total        = (2c) / (2t) = c/t  ✓
        #   off-diagonal: counts[i,j] / (2 * total)  = (2c) / (4t) = c/(2t) ✓
        #
        # Without the ÷2 for off-diagonal, those entries would be 2× too large.
        mat = zeros(Rational{Int}, nf, nf)
        for i in 1:nf, j in i:nf
            mat[i, j] = if total == 0
                0//1
            elseif i == j
                counts[i, j] // total
            else
                counts[i, j] // (2 * total)
            end
        end
        mat
    end
end

# Given a fixed type injection and a set of unlabeled vertices,
# find the index of the corresponding flag in flag_canons.
function find_flag(H::Hypergraph{R}, type_perm, half_verts,
                   flag_canons::Vector{Vector{NTuple{R,Int}}},
                   s::Int) :: Union{Int,Nothing} where R
    verts     = [collect(type_perm); collect(half_verts)]
    sub       = induce(H, verts)
    sub_canon = canonical(Flag{R}(sub, s)).graph.edges
    findfirst(c -> c == sub_canon, flag_canons)
end
