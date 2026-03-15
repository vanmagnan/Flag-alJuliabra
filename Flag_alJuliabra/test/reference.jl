# Utilities for parsing C flagmatic reference output.
# Used only in tests — not part of the FlagAlgebra package.

using FlagAlgebra

# --- C graph string parser --------------------------------------------------

# Parse the C graph string format into a Hypergraph.
# Format: "n:e1e2...em" where each edge is K consecutive digit characters.
# Example (K=2): "4:121314"  → 4 vertices, edges (1,2),(1,3),(1,4)
# Example (K=3): "5:123234"  → 5 vertices, 3-edges (1,2,3),(2,3,4)
# Special case:  "0:"        → empty hypergraph
function parse_c_graph(s::AbstractString, K::Int; oriented::Bool=false) :: Hypergraph
    sep = findfirst(':', s)
    sep === nothing && error("Invalid graph string: $s")
    n = parse(Int, s[1:sep-1])
    edge_str = s[sep+1:end]
    isempty(edge_str) && return Hypergraph{K}(n, oriented, NTuple{K,Int}[])
    length(edge_str) % K == 0 || error("Edge string length not divisible by K=$K: $s")
    edges = [NTuple{K,Int}(parse(Int, string(edge_str[i+j])) for j in 0:K-1)
             for i in 1:K:length(edge_str)]
    Hypergraph{K}(n, oriented, edges)
end

# Parse a graph string into a Flag with a specified type_size.
function parse_c_flag(s::AbstractString, K::Int, type_size::Int; oriented::Bool=false) :: Flag
    Flag{K}(parse_c_graph(s, K; oriented), type_size)
end

# --- flags.py parser --------------------------------------------------------

# Extract a named list from flags.py, e.g. types = ["...", "...", ...]
# Returns the raw string content between the outermost brackets.
function extract_list(text::AbstractString, name::AbstractString) :: Vector{String}
    # Match: name = [\n    "...",\n    "..."\n]
    pattern = Regex("$name\\s*=\\s*\\[([^\\[\\]]*?)\\]", "s")
    m = match(pattern, text)
    m === nothing && error("Could not find '$name' in flags.py")
    # Extract quoted strings
    [String(m2.match[2:end-1]) for m2 in eachmatch(r"\"[^\"]*\"", m.captures[1])]
end

# Extract a scalar integer from flags.py.
function extract_int(text::AbstractString, name::AbstractString) :: Int
    m = match(Regex("$name\\s*=\\s*(\\d+)"), text)
    m === nothing && error("Could not find '$name' in flags.py")
    parse(Int, m.captures[1])
end

# Extract num_flags list from flags.py (it's a list of ints, not strings).
function extract_int_list(text::AbstractString, name::AbstractString) :: Vector{Int}
    m = match(Regex("$name\\s*=\\s*\\[([^\\]]*)\\]"), text)
    m === nothing && error("Could not find '$name' in flags.py")
    [parse(Int, strip(s)) for s in split(m.captures[1], ",") if !isempty(strip(s))]
end

# Extract the nested flags list: flags = [ [...], [...], ... ]
# Returns a Vector of Vectors of graph strings.
# Uses a line-by-line approach to avoid fragile byte-offset arithmetic.
function extract_nested_list(text::AbstractString) :: Vector{Vector{String}}
    result   = Vector{Vector{String}}()
    current  = String[]
    in_flags = false   # have we entered the outer "flags = [" block yet
    in_inner = false   # are we currently inside an inner [...] block

    for line in split(text, '\n')
        s = strip(line)
        if !in_flags
            # Wait until we see the "flags = [" line
            (s == "flags = [" || startswith(s, "flags = [")) && (in_flags = true)
        elseif s == "[" || s == "],["   # opening of an inner block
            in_inner = true
            current = String[]
        elseif (s == "]" || s == "],") && in_inner
            # closing of an inner block — save it and reset
            push!(result, current)
            current = String[]
            in_inner = false
        elseif s == "]"                 # closing of the outer flags block — done
            break
        elseif in_inner
            # inside an inner block: extract any quoted graph string on this line
            m = match(r"\"([^\"]*)\"", s)
            m !== nothing && push!(current, String(m.captures[1]))
        end
    end
    result
end

# Parse a complete flags.py file into a named tuple.
function parse_flags_py(path::AbstractString, K::Int) :: NamedTuple
    text = read(path, String)
    n          = extract_int(text, "n")
    num_types  = extract_int(text, "num_types")
    num_flags  = extract_int_list(text, "num_flags")
    type_strs  = extract_list(text, "types")
    flag_strs  = extract_nested_list(text)
    H_strs     = extract_list(text, "H")

    # Type size is the vertex count of the type graph
    types = [parse_c_flag(s, K, parse_c_graph(s, K).n) for s in type_strs]
    flags = [[parse_c_flag(s, K, types[σ].type_size) for s in flag_strs[σ]]
             for σ in 1:num_types]
    H     = [parse_c_graph(s, K) for s in H_strs]

    (n=n, num_types=num_types, num_flags=num_flags, types=types, flags=flags, H=H)
end

# --- flags.rat parser -------------------------------------------------------

# Parse flags.rat into a lookup dictionary.
# Format per line: H_idx block_idx flag_i flag_j numer denom
# block_idx = σ+1 (type block 2 = type 1, etc.)
#
# Returns Dict mapping (H_idx, σ, i, j) => Rational{Int}
# where σ is 1-based type index, i ≤ j are 1-based flag indices.
function parse_flags_rat(path::AbstractString) :: Dict{NTuple{4,Int}, Rational{Int}}
    result = Dict{NTuple{4,Int}, Rational{Int}}()
    for line in eachline(path)
        parts = split(strip(line))
        length(parts) == 6 || continue
        H_idx, block, i, j, numer, denom = parse.(Int, parts)
        σ = block - 1   # block 2 = type 1
        result[(H_idx, σ, i, j)] = numer // denom
    end
    result
end
