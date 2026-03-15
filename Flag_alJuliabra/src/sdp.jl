# SDP assembly via JuMP.
# The solver is not specified here — the user attaches one after calling build_sdp.
#
# Example:
#   using COSMO
#   model = build_sdp(data)
#   set_optimizer(model, COSMO.Optimizer)
#   optimize!(model)
#   objective_value(model)

using JuMP
using LinearAlgebra

# Build the flag algebra SDP for a FlagAlgebraData.
#
# Two modes, controlled by prob.minimize:
#
#   Upper bound (minimize=false, maximize density):
#     Minimize λ  subject to:  λ - Σ_σ <P_σ(H_i), Q_σ> - s_i = density(H_i)
#     Since Q_σ ≽ 0 and P_σ(H_i) ≽ 0 → flag_sum ≥ 0, slack ≥ 0 → λ ≥ density(H_i).
#     Minimizing λ gives the tightest upper bound.
#
#   Lower bound (minimize=true, minimize density):
#     Maximize λ  subject to:  λ + Σ_σ <P_σ(H_i), Q_σ> - s_i = density(H_i)
#     flag_sum ≥ 0, slack ≥ 0 → λ ≤ density(H_i).
#     Maximizing λ gives the tightest lower bound.
#
# In both cases the constraint is:
#   λ + sign * flag_sum - s_i = density(H_i)   (sign = +1 for lower bound, -1 for upper)
#
# Sharp graphs (i ∈ sharps) have s_i = 0 (equality, no slack variable).
#
# The inner product <P, Q> for symmetric Q and upper-triangular P expands to:
#   Σ_j P[j,j]*Q[j,j]  +  2 * Σ_{j<k} P[j,k]*Q[j,k]

# build_sdp returns a named tuple (model, Q) so that the Q variables remain
# accessible after solving for certificate extraction.
function build_sdp(data::FlagAlgebraData{R};
                   sharps::Vector{Int} = Int[]) where R

    model = Model()

    @variable(model, λ)

    # Upper bound: minimize λ (λ ≥ densities).
    # Lower bound: maximize λ (λ ≤ densities).
    if data.problem.minimize
        @objective(model, Max, λ)
    else
        @objective(model, Min, λ)
    end

    # One PSD matrix per type
    Q = map(data.flags) do type_flags
        nf = length(type_flags)
        @variable(model, [1:nf, 1:nf], PSD)
    end

    # Slack variables for non-sharp admissible graphs
    num_H = length(data.admissible)
    @variable(model, slack[1:num_H] >= 0)

    # Constraint sign: -1 for upper bound (flag_sum subtracted), +1 for lower bound (added).
    csign = data.problem.minimize ? 1 : -1

    for (i, _) in enumerate(data.admissible)
        # Build Σ_σ <P_σ(H_i), Q_σ> as an affine JuMP expression.
        # Inner product for upper-triangular P and symmetric Q:
        #   diagonal (j==k): P[j,j] * Q[j,j]
        #   off-diagonal (j<k): P[j,k] * 2 * Q[j,k]  (factor-of-2 from symmetry of Q)
        flag_sum = AffExpr(0.0)
        for σ in eachindex(data.types)
            nf = length(data.flags[σ])
            P  = data.pair_dens[i][σ]
            for j in 1:nf, k in j:nf
                P[j, k] == 0 && continue
                add_to_expression!(flag_sum,
                    Float64(P[j, k]) * (j == k ? 1 : 2),
                    Q[σ][j, k])
            end
        end

        if i ∈ sharps
            @constraint(model, λ + csign * flag_sum == Float64(data.densities[i]))
        else
            @constraint(model, λ + csign * flag_sum - slack[i] == Float64(data.densities[i]))
        end
    end

    (model=model, Q=Q)
end

# Attach a solver, optimize, and return a FlagAlgebraResult.
#
# optimizer   — a JuMP-compatible optimizer constructor, e.g. COSMO.Optimizer
# extract_Q   — if true, populate result.Q with the certificate PSD matrices
#
# The solver is silenced by default; pass silent=false to see solver output.
function solve_sdp(data::FlagAlgebraData{R}, optimizer;
                   sharps::Vector{Int}  = Int[],
                   extract_Q::Bool      = false,
                   silent::Bool         = true) :: FlagAlgebraResult{R} where R

    built = build_sdp(data; sharps=sharps)
    set_optimizer(built.model, optimizer)
    silent && set_silent(built.model)
    optimize!(built.model)

    status = Symbol(termination_status(built.model))
    bound  = objective_value(built.model)
    Q_vals = extract_Q ? [value.(q) for q in built.Q] : nothing

    FlagAlgebraResult{R}(data.problem, status, bound, Q_vals)
end

# Verify a flag algebra certificate in exact arithmetic.
#
# Steps:
#   1. Rationalize λ and each Q matrix (Float64 → Rational{Int64}).
#   2. Check each Q is numerically PSD (min eigenvalue ≥ -psd_tol).
#   3. For each admissible graph H_i, compute the exact rational residual:
#        upper bound:  residual[i] = λ_rat - Σ_σ <Q_σ, P_σ(H_i)> - density[i]
#        lower bound:  residual[i] = density[i] - λ_rat - Σ_σ <Q_σ, P_σ(H_i)>
#      A non-negative residual certifies the bound for H_i.
#
# Returns a named tuple:
#   valid          — true iff all Q PSD (within psd_tol) and all residuals ≥ 0
#   λ_certified    — rationalized bound
#   min_psd_eigval — minimum eigenvalue across all Q matrices (Float64)
#   residuals      — rational residual per admissible graph
#   min_residual   — tightest (smallest) residual
#
# rat_tol controls the rationalize approximation: smaller = more exact but
# larger denominators. For typical solver accuracy (~1e-6), rat_tol=1e-6 works.
# Overflow is possible for large problems; use rat_tol ≥ 1e-8 to stay safe.

function verify_certificate(data::FlagAlgebraData{R}, result::FlagAlgebraResult{R};
                             rat_tol::Float64 = 1e-6,
                             psd_tol::Float64 = 1e-6) where R

    result.Q === nothing &&
        error("No Q matrices in result — rerun solve_sdp with extract_Q=true")

    minimize = data.problem.minimize

    # Step 1: rationalize λ and Q matrices.
    # Q uses BigInt to avoid overflow when pair-density and Q entries are accumulated —
    # the rational LCM of many small fractions can exceed Int64/Int128.
    λ_rat = rationalize(Int, result.bound; tol=rat_tol)
    Q_rat = [rationalize.(BigInt, q; tol=rat_tol) for q in result.Q]

    # Step 2: PSD check — use Float64 eigenvalues on each Q (Symmetric for stability)
    min_psd_eigval = minimum(
        minimum(eigvals(Symmetric(Float64.(q))))
        for q in result.Q
    )
    psd_ok = min_psd_eigval ≥ -psd_tol

    # Step 3: exact rational residuals
    # Inner product <Q_σ, P_σ(H_i)>: diagonal entries contribute once,
    # off-diagonal entries contribute twice (factor-of-2 from symmetry of Q).
    residuals = map(enumerate(data.admissible)) do (i, _)
        flag_sum = zero(Rational{BigInt})
        for σ in eachindex(data.types)
            nf  = length(data.flags[σ])
            P   = data.pair_dens[i][σ]
            Q_s = Q_rat[σ]
            for j in 1:nf, k in j:nf
                P[j, k] == 0 && continue
                flag_sum += P[j, k] * (j == k ? 1 : 2) * Q_s[j, k]
            end
        end

        minimize ?
            data.densities[i] - λ_rat - flag_sum :   # lower bound: density - λ - sum ≥ 0
            λ_rat - flag_sum - data.densities[i]      # upper bound: λ - sum - density ≥ 0
    end

    min_residual = minimum(residuals)
    valid = psd_ok && min_residual ≥ 0

    (valid          = valid,
     λ_certified    = λ_rat,
     min_psd_eigval = min_psd_eigval,
     residuals      = residuals,
     min_residual   = min_residual)
end

# Identify sharp (extremal) admissible graphs from a certificate.
#
# Sharp graphs have zero slack: residual == 0 exactly (as a rational).
# These are the graphs where the bound is tight — the extremal examples.
#
# Returns a named tuple:
#   indices   — positions in data.admissible with residual == 0
#   graphs    — the corresponding Hypergraph objects
#   densities — their densities
#   residuals — all residuals (for inspection)
function identify_sharps(data::FlagAlgebraData{R}, result::FlagAlgebraResult{R};
                         rat_tol::Float64 = 1e-6,
                         psd_tol::Float64 = 1e-6) where R
    cert = verify_certificate(data, result; rat_tol=rat_tol, psd_tol=psd_tol)
    indices = findall(r -> r == 0//1, cert.residuals)
    (indices   = indices,
     graphs    = data.admissible[indices],
     densities = data.densities[indices],
     residuals = cert.residuals)
end
