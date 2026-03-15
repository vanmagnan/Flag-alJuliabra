# Flag and admissible graph generation.

using Combinatorics

# All possible R-edges on vertices 1:n (sorted tuples).
all_possible_edges(n::Int, ::Val{R}) where R =
    NTuple{R,Int}.(combinations(1:n, R))

# Generate all R-uniform flags on n vertices over a given type,
# subject to forbidden subgraph constraints.
#
# A flag is a hypergraph where the first type.type_size vertices are labeled.
# The type itself is the seed: we add edges one layer at a time,
# only extending graphs from the previous layer (those with e-1 extra edges).
#
# forbidden         : no subgraph of this kind is allowed
# forbidden_induced : no induced subgraph of this kind is allowed (checked post-generation)

function generate_flags(n::Int, type::Flag{R},
                        forbidden::Vector{Hypergraph{R}},
                        forbidden_induced::Vector{Hypergraph{R}}) :: Vector{Flag{R}} where R

    s = type.type_size

    # Edges that involve at least one unlabeled vertex.
    # Edges entirely within the type are already present and cannot be added.
    candidate_edges = filter(all_possible_edges(n, Val(R))) do e
        any(v -> v > s, e)
    end

    # Seed flag: type edges embedded in an n-vertex graph, so the 0-extra-edges
    # flag (type extended with all unlabeled vertices isolated) is included in output.
    seed_graph = Hypergraph{R}(n, type.graph.oriented, type.graph.edges)
    seed_flag  = Flag{R}(seed_graph, s)

    # layers[i] = flags with exactly i-1 extra edges beyond the type.
    # Using Vector for cache-friendly iteration (mirrors C's e_start[] approach).
    layers = Vector{Vector{Flag{R}}}([[seed_flag]])

    for _ in 1:length(candidate_edges)
        current_layer   = Vector{Flag{R}}()
        seen_this_layer = Set{Vector{NTuple{R,Int}}}()  # keyed on canonical edge list

        for f in layers[end]
            edge_set = Set(f.graph.edges)

            for new_edge in candidate_edges
                new_edge ∈ edge_set && continue

                new_edges = sort([f.graph.edges; [new_edge]])
                g         = Hypergraph{R}(n, f.graph.oriented, new_edges)
                candidate = Flag{R}(g, s)

                # Check forbidden subgraphs only on vertex subsets containing new_edge.
                # This mirrors the C optimization: only check combinations that include
                # the newly added edge, since all smaller subsets were clean before.
                is_bad = any(forbidden) do fg
                    any(combinations(1:n, fg.n)) do verts
                        !issubset(new_edge, verts) && return false
                        sub = induce(g, collect(verts))
                        has_subgraph(sub, fg)
                    end
                end
                is_bad && continue

                c = canonical(candidate).graph.edges
                if c ∉ seen_this_layer
                    push!(seen_this_layer, c)
                    push!(current_layer, canonical(candidate))
                end
            end
        end

        isempty(current_layer) && break
        push!(layers, current_layer)
    end

    # Include layers[1] (the 0-extra-edges seed flag) through all generated layers.
    all_flags = reduce(vcat, layers; init = Flag{R}[])

    # Second pass: remove flags containing forbidden induced subgraphs.
    filter!(all_flags) do f
        !any(forbidden_induced) do fg
            has_induced_subgraph(f.graph, fg)
        end
    end
end

# Generate admissible graphs: flags with no labeled vertices (type_size = 0).
function generate_admissible(n::Int, ::Val{R};
                              forbidden::Vector{Hypergraph{R}}         = Hypergraph{R}[],
                              forbidden_induced::Vector{Hypergraph{R}} = Hypergraph{R}[],
                              oriented::Bool                            = false) :: Vector{Hypergraph{R}} where R
    seed = Flag{R}(empty_hypergraph(R, oriented), 0)
    # Admissible graphs are flags over the empty type; re-wrap as plain Hypergraph
    flags = generate_flags(n, seed, forbidden, forbidden_induced)
    map(f -> f.graph, flags)
end

# Generate types of a given order: admissible graphs with all vertices labeled.
function generate_types(order::Int, ::Val{R};
                        forbidden::Vector{Hypergraph{R}}         = Hypergraph{R}[],
                        forbidden_induced::Vector{Hypergraph{R}} = Hypergraph{R}[],
                        oriented::Bool                            = false) :: Vector{Flag{R}} where R
    graphs = generate_admissible(order, Val(R);
                                 forbidden = forbidden,
                                 forbidden_induced = forbidden_induced,
                                 oriented = oriented)
    # A type is a flag where all vertices are labeled
    [Flag{R}(g, g.n) for g in graphs]
end
