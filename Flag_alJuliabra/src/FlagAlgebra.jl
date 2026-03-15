module FlagAlgebra

include("types.jl")
include("graphs.jl")
include("generation.jl")
include("densities.jl")
include("sdp.jl")
include("pipeline.jl")

export Hypergraph, Flag, FlagProblem, FlagAlgebraData, FlagAlgebraResult
export Graph, Threegraph, Fourgraph
export complete_hypergraph, k4_minus, c5_3uniform, f32
export induce, canonical, flag_isomorphic, has_subgraph, has_induced_subgraph
export generate_flags, generate_admissible, generate_types
export edge_density, induced_density, compute_pair_densities
export build_sdp, solve_sdp, verify_certificate, identify_sharps
export build_flag_algebra_data

end
