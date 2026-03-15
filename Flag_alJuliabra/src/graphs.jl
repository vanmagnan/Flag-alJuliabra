# Graph operations: induction, canonical forms, subgraph tests.

using Combinatorics

# --- Induced subgraph -------------------------------------------------------

# Induce the subgraph of g on the given vertex subset (1-based indices).
# Vertices are reindexed 1:length(verts) in the order given.
function induce(g::Hypergraph{R}, verts::AbstractVector{Int}) :: Hypergraph{R} where R
    idx = Dict(v => i for (i, v) in enumerate(verts))
    new_edges = NTuple{R,Int}[]
    for e in g.edges
        if all(v -> haskey(idx, v), e)
            push!(new_edges, normalize_edge(NTuple{R,Int}(idx[v] for v in e)))
        end
    end
    Hypergraph{R}(length(verts), g.oriented, sort(new_edges))
end

# --- Canonical form ---------------------------------------------------------

# Compute the canonical form of a flag by minimizing the edge set
# over all permutations of the unlabeled vertices (indices type_size+1 : n).
# Labeled vertices (1 : type_size) are held fixed.
function canonical(f::Flag{R}) :: Flag{R} where R
    s, n = f.type_size, f.graph.n
    best = f

    for perm in permutations(s+1:n)
        relabel = Dict{Int,Int}(i => i for i in 1:s)
        for (j, v) in enumerate(perm)
            relabel[s + j] = v
        end
        new_edges = sort([normalize_edge(NTuple{R,Int}(relabel[v] for v in e))
                          for e in f.graph.edges])
        candidate = Flag{R}(Hypergraph{R}(n, f.graph.oriented, new_edges), s)
        if new_edges < best.graph.edges
            best = candidate
        end
    end
    best
end

# Two flags are isomorphic (as flags, fixing labeled vertices) iff their
# canonical forms are equal.
function flag_isomorphic(f1::Flag{R}, f2::Flag{R}) :: Bool where R
    f1.type_size == f2.type_size &&
    f1.graph.n   == f2.graph.n   &&
    canonical(f1).graph.edges == canonical(f2).graph.edges
end

# --- Subgraph tests ---------------------------------------------------------

# Does g contain sg as a (not necessarily induced) subgraph?
# Checks all injections of V(sg) into V(g).
function has_subgraph(g::Hypergraph{R}, sg::Hypergraph{R}) :: Bool where R
    edge_set = Set(g.edges)
    for perm in permutations(1:g.n, sg.n)
        if all(normalize_edge(NTuple{R,Int}(perm[v] for v in e)) ∈ edge_set
               for e in sg.edges)
            return true
        end
    end
    false
end

# Does g contain sg as an induced subgraph?
function has_induced_subgraph(g::Hypergraph{R}, sg::Hypergraph{R}) :: Bool where R
    sg_canon = canonical(Flag{R}(sg, 0))
    for verts in combinations(1:g.n, sg.n)
        sub = induce(g, verts)
        canonical(Flag{R}(sub, 0)).graph.edges == sg_canon.graph.edges && return true
    end
    false
end
