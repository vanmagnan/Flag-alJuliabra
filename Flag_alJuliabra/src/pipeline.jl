# High-level pipeline: FlagProblem → FlagAlgebraData.
# Ties together generation, density computation, and struct assembly.

# Convenience constructor for FlagProblem — keyword arguments for optional fields.
function FlagProblem(n::Int, type_order::Int, ::Val{R};
                     forbidden::Vector{Hypergraph{R}}         = Hypergraph{R}[],
                     forbidden_induced::Vector{Hypergraph{R}} = Hypergraph{R}[],
                     target::Union{Hypergraph{R}, Nothing}    = nothing,
                     minimize::Bool                           = false) where R
    FlagProblem{R}(n, type_order, forbidden, forbidden_induced, target, minimize)
end

# Run the full pre-SDP pipeline for a FlagProblem.
#
# Steps:
#   1. generate_types  — all non-isomorphic types at each valid order
#   2. generate_flags  — flags over each type at the appropriate flag size
#   3. generate_admissible — admissible graphs on n vertices
#   4. densities       — edge density (or induced density of target) for each H
#   5. compute_pair_densities — pair density matrices for every (H, type) pair
#
# Valid type orders: all s with n-s even (same parity as n),
#   from n%2 up to prob.type_order in steps of 2.

function build_flag_algebra_data(prob::FlagProblem{R}) :: FlagAlgebraData{R} where R
    n = prob.n

    # Collect types from all valid orders
    min_s = n % 2   # 0 if n even, 1 if n odd (so that n-s is always even)
    types = reduce(vcat, [
        generate_types(s, Val(R);
                       forbidden         = prob.forbidden,
                       forbidden_induced = prob.forbidden_induced)
        for s in min_s:2:prob.type_order
    ])

    # Flags over each type at flag size m = ⌊(n+s)/2⌋
    flags = map(types) do t
        s = t.type_size
        m = isodd(n - s) ? (n + s - 1) ÷ 2 : (n + s) ÷ 2
        generate_flags(m, t, prob.forbidden, prob.forbidden_induced)
    end

    # Admissible graphs on n vertices
    admissible = generate_admissible(n, Val(R);
                                     forbidden         = prob.forbidden,
                                     forbidden_induced = prob.forbidden_induced)

    # Density of each admissible graph
    densities = map(admissible) do H
        prob.target === nothing ? edge_density(H) : induced_density(H, prob.target)
    end

    # Pair densities: pair_dens[H_idx][σ] = upper-triangular matrix for that (H, type) block
    # Parallelized over admissible graphs — each call is independent (no shared mutable state).
    pair_dens = Vector{Vector{Matrix{Rational{Int}}}}(undef, length(admissible))
    Threads.@threads :static for i in eachindex(admissible)
        pair_dens[i] = compute_pair_densities(admissible[i], types, flags)
    end

    FlagAlgebraData{R}(prob, types, flags, admissible, densities, pair_dens)
end
