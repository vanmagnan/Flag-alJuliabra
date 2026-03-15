# Core data structures for the flag algebra method.
#
# An R-uniform hypergraph has every edge of size R.
# R=2 gives ordinary graphs; R=3 gives 3-uniform hypergraphs.

using Combinatorics

# Edges are stored as sorted R-tuples of 1-based vertex indices.
# Sorting is the canonical form for a single edge (undirected, unoriented).

struct Hypergraph{R}
    n        :: Int                      # number of vertices
    oriented :: Bool
    edges    :: Vector{NTuple{R,Int}}    # sorted R-tuples; sorted as a set too

    # Inner constructor normalizes and sorts edges, then calls new{R} directly
    # to avoid recursion. All construction goes through here.
    function Hypergraph{R}(n::Int, oriented::Bool, edges::AbstractVector) where R
        sorted = sort(normalize_edge.(NTuple{R,Int}.(edges)))
        new{R}(n, oriented, sorted)
    end
end

# A Flag is a hypergraph where the first type_size vertices are labeled.
# The labeled vertices are fixed when computing canonical forms —
# only the unlabeled vertices (type_size+1 : n) are permuted.
#
# Special cases:
#   type_size == 0         → admissible graph (no labeled vertices)
#   type_size == graph.n   → type (all vertices labeled)

struct Flag{R}
    graph     :: Hypergraph{R}
    type_size :: Int
end

# A FlagProblem specifies which SDP to build: the uniformity, the vertex
# count for admissible graphs, and any constraints on which graphs appear.

struct FlagProblem{R}
    n                 :: Int                       # vertex count for admissible graphs
    type_order        :: Int                       # vertex count for types (must be ≤ n-2, same parity as n)
    forbidden         :: Vector{Hypergraph{R}}     # forbidden subgraphs
    forbidden_induced :: Vector{Hypergraph{R}}     # forbidden induced subgraphs
    target            :: Union{Hypergraph{R}, Nothing}  # graph to optimize density of; nothing = edge density
    minimize          :: Bool
end

# FlagAlgebraData holds everything computed from a FlagProblem.
# This is what gets passed to build_sdp.

struct FlagAlgebraData{R}
    problem    :: FlagProblem{R}
    types      :: Vector{Flag{R}}
    flags      :: Vector{Vector{Flag{R}}}                    # flags[σ] = flags over types[σ]
    admissible :: Vector{Hypergraph{R}}
    densities  :: Vector{Rational{Int}}                      # density of each admissible graph
    pair_dens  :: Vector{Vector{Matrix{Rational{Int}}}}      # pair_dens[H_idx][σ][i,j] (upper triangular)
end

# FlagAlgebraResult holds the output of solve_sdp.

struct FlagAlgebraResult{R}
    problem :: FlagProblem{R}
    status  :: Symbol                                   # JuMP termination status as a Symbol
    bound   :: Float64                                  # optimal λ
    Q       :: Union{Vector{Matrix{Float64}}, Nothing}  # PSD certificate matrices (one per type), or nothing
end

# --- Constructors -----------------------------------------------------------

# Normalize a single edge: sort vertices (for undirected hypergraphs)
normalize_edge(e::NTuple{R,Int}) where R = NTuple{R,Int}(sort(collect(e)))

# Convenience: unoriented graph/hypergraph
Graph(n::Int, edges)      = Hypergraph{2}(n, false, edges)
Threegraph(n::Int, edges) = Hypergraph{3}(n, false, edges)
Fourgraph(n::Int, edges)  = Hypergraph{4}(n, false, edges)

# Convenience: empty hypergraph (used as seed type in generation)
empty_hypergraph(r::Int, oriented::Bool=false) = Hypergraph{r}(0, oriented, NTuple{r,Int}[])

# --- Named graph shortcuts -----------------------------------------------
#
# Mirrors the C flagmatic CLI shortcuts (--forbid-k4-, --forbid-k4, etc.)
# so users don't need to write out edge lists for common forbidden graphs.

# Complete R-uniform hypergraph on n vertices: all C(n,R) edges.
# Covers --forbid-k4  (complete_hypergraph(4, Val(3)))
#    and --forbid-k5  (complete_hypergraph(5, Val(3)))
#    and R=4 analogs  (complete_hypergraph(5, Val(4)), etc.)
complete_hypergraph(n::Int, ::Val{R}) where R =
    Hypergraph{R}(n, false, NTuple{R,Int}.(combinations(1:n, R)))

# K₄⁻: the unique 3-graph on 4 vertices with 3 edges (K₄ minus one edge).
# Corresponds to C's --forbid-k4-  ("4.3").
k4_minus() = Threegraph(4, [(1,2,3),(1,2,4),(1,3,4)])

# Tight (linear) 5-cycle in 3-uniform hypergraphs.
# Edges: {1,2,3},{2,3,4},{3,4,5},{4,5,1},{5,1,2}.
# Corresponds to C's --forbid-c5  ("5:123234345451512").
c5_3uniform() = Threegraph(5, [(1,2,3),(2,3,4),(3,4,5),(4,5,1),(5,1,2)])

# F₃₂: 5-vertex 3-graph with 4 edges sharing a common pair {4,5}.
# Edges: {1,2,3},{1,4,5},{2,4,5},{3,4,5}.
# Corresponds to C's --forbid-f32  ("5:123145245345").
f32() = Threegraph(5, [(1,2,3),(1,4,5),(2,4,5),(3,4,5)])

# --- Display ----------------------------------------------------------------

function Base.show(io::IO, g::Hypergraph{R}) where R
    print(io, "$(g.n)v $(R)-uniform hypergraph, $(length(g.edges)) edges")
end

function Base.show(io::IO, res::FlagAlgebraResult{R}) where R
    println(io, "FlagAlgebraResult{$R}")
    println(io, "  status : $(res.status)")
    println(io, "  bound  : $(res.bound)")
    print(  io, "  Q      : ", res.Q === nothing ? "not extracted" : "$(length(res.Q)) matrices")
end

function Base.show(io::IO, f::Flag{R}) where R
    if f.type_size == 0
        print(io, "Admissible: $(f.graph)")
    elseif f.type_size == f.graph.n
        print(io, "Type: $(f.graph)")
    else
        print(io, "Flag over $(f.type_size)-vertex type: $(f.graph)")
    end
end
