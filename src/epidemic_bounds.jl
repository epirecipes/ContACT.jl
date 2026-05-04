"""
Epidemic bounds from partial next-generation matrix information.

Implements results from Britton, Poletti, Scarpa & Pellis (2025):
"Bounds on R₀ and final epidemic size when the next-generation matrix M
is only partially known" (arXiv:2602.23885v2).

Convention: The NGM here follows the paper's convention where M[i,j] is the
expected number of infectious contacts that an i-individual has with
j-individuals. Row sums rᵢ = Σⱼ M[i,j] = total infectious output of type i.
This is the TRANSPOSE of ContACT's internal next_generation_matrix (which is
K[target, source]).
"""

# ---------------------------------------------------------------------------
# Scalar final-size solvers
# ---------------------------------------------------------------------------

"""
    solve_final_size_scalar(α)

Compute tₐ, the largest solution to `1 - t = exp(-α·t)` in [0,1].

Returns 0 when α ≤ 1 (subcritical). For α > 1, uses bisection on the
unique positive root.
"""
function solve_final_size_scalar(α::Real)
    α = Float64(α)
    isfinite(α) && α >= 0 || throw(ArgumentError("α must be finite and non-negative"))
    α <= 1.0 && return 0.0

    # f(t) = 1 - t - exp(-αt); f(0)=0, f(ε)>0 for small ε when α>1, f(1)<0
    # The positive root is in (0, 1).
    lo, hi = 1e-10, 1.0 - 1e-15
    for _ in 1:100
        mid = (lo + hi) / 2
        val = 1.0 - mid - exp(-α * mid)
        if val > 0
            lo = mid
        else
            hi = mid
        end
        (hi - lo) < 1e-14 && break
    end
    (lo + hi) / 2
end

"""
    solve_final_size_ext(α, γ)

Compute t_{α,γ}, the largest solution to `t = 1 - exp(-α·t - γ)` in [0,1].

This extends the basic final-size equation with an external force of infection γ.
For γ = 0, this reduces to `solve_final_size_scalar(α)`.
"""
function solve_final_size_ext(α::Real, γ::Real)
    α = Float64(α)
    γ = Float64(γ)
    isfinite(α) && α >= 0 || throw(ArgumentError("α must be finite and non-negative"))
    isfinite(γ) && γ >= 0 || throw(ArgumentError("γ must be finite and non-negative"))

    # g(t) = 1 - exp(-αt - γ) - t; find largest root in [0,1]
    # g(0) = 1 - exp(-γ) ≥ 0 (with equality iff γ=0)
    # g(1) = 1 - exp(-α - γ) - 1 = -exp(-α - γ) < 0
    # When γ > 0, there is always a positive root.
    # When γ = 0 and α ≤ 1, only root is 0.
    if γ == 0.0
        return solve_final_size_scalar(α)
    end

    # With γ > 0, g(0) > 0 and g(1) < 0, so root in (0,1)
    lo, hi = 0.0, 1.0
    for _ in 1:100
        mid = (lo + hi) / 2
        val = 1.0 - exp(-α * mid - γ) - mid
        if val > 0
            lo = mid
        else
            hi = mid
        end
        (hi - lo) < 1e-14 && break
    end
    (lo + hi) / 2
end

"""
    solve_final_size_vector(K, π)

Solve the multitype final-size equation `1 - τⱼ = exp(-Σᵢ πᵢτᵢ Kᵢⱼ/πⱼ)`
for the vector τ (paper convention: K[i,j] = infections from i to j-individuals).

Returns the largest non-trivial solution via fixed-point iteration.
Returns zeros if R₀ ≤ 1.
"""
function solve_final_size_vector(K::AbstractMatrix{<:Real}, π::AbstractVector{<:Real})
    k = length(π)
    size(K) == (k, k) || throw(DimensionMismatch("K must be $k × $k"))
    all(x -> isfinite(x) && x >= 0, K) || throw(ArgumentError("K must be non-negative"))
    all(x -> isfinite(x) && x >= 0, π) || throw(ArgumentError("π must be non-negative"))

    # Check R₀
    R0 = maximum(abs.(eigvals(K)))
    R0 > 1.0 || return zeros(k)

    # Fixed-point iteration: τⱼ ← 1 - exp(-Σᵢ πᵢτᵢ Kᵢⱼ/πⱼ)
    # Start near the scalar solution
    τ = fill(solve_final_size_scalar(R0), k)
    for _ in 1:10000
        τ_new = zeros(k)
        for j in 1:k
            if π[j] == 0
                τ_new[j] = 0.0
                continue
            end
            exponent = sum(π[i] * τ[i] * K[i, j] / π[j] for i in 1:k if π[i] > 0)
            τ_new[j] = 1.0 - exp(-exponent)
        end
        if maximum(abs.(τ_new .- τ)) < 1e-12
            return τ_new
        end
        τ = τ_new
    end
    τ
end

# ---------------------------------------------------------------------------
# R₀ bounds
# ---------------------------------------------------------------------------

"""
    EpidemicBounds

Result type for epidemic bounds, holding lower and upper values.
"""
struct EpidemicBounds{T}
    lower::T
    upper::T
end

Base.show(io::IO, b::EpidemicBounds) = print(io, "[$(b.lower), $(b.upper)]")

"""
    r0_bounds(K; info=:row)

Compute sharp R₀ bounds for a general next-generation matrix when only row
or column sums are known.

# Arguments
- `K`: Next-generation matrix (paper convention: K[i,j] = infections from i to j)
- `info`: `:row` (row sums known), `:col` (column sums known), or `:both`

# Returns
`EpidemicBounds` with lower and upper R₀ bounds.

Based on Theorem 3.1 of Britton et al. (2025).
"""
function r0_bounds(K::AbstractMatrix{<:Real}; info::Symbol=:row)
    info in (:row, :col, :both) || throw(ArgumentError("info must be :row, :col, or :both"))
    k = size(K, 1)
    size(K, 2) == k || throw(DimensionMismatch("K must be square"))
    all(x -> isfinite(x) && x >= 0, K) || throw(ArgumentError("K must be non-negative"))

    row_sums = vec(sum(K; dims=2))  # rᵢ = Σⱼ K[i,j]
    col_sums = vec(sum(K; dims=1))  # cⱼ = Σᵢ K[i,j]

    if info == :row
        return EpidemicBounds(minimum(row_sums), maximum(row_sums))
    elseif info == :col
        return EpidemicBounds(minimum(col_sums), maximum(col_sums))
    else  # :both
        lo = max(minimum(row_sums), minimum(col_sums))
        hi = min(maximum(row_sums), maximum(col_sums))
        return EpidemicBounds(lo, hi)
    end
end

"""
    r0_bounds(cm::ContactMatrix; info=:row, transmissibility=1, recovery_rate=1)

Convenience wrapper that computes the NGM from a `ContactMatrix` and returns
R₀ bounds. The NGM is transposed to paper convention internally.
"""
function r0_bounds(cm::ContactMatrix; info::Symbol=:row,
                   transmissibility::Real=1, recovery_rate::Real=1)
    K_contact = next_generation_matrix(cm;
        transmissibility=transmissibility, recovery_rate=recovery_rate)
    K_paper = transpose(K_contact)
    r0_bounds(Matrix(K_paper); info=info)
end

"""
    r0_bounds_detailed_balance(K, π; info=:row)

Compute R₀ bounds when the NGM satisfies detailed balance (πᵢKᵢⱼ = πⱼKⱼᵢ).

# Arguments
- `K`: Next-generation matrix (paper convention)
- `π`: Population fractions (must sum to ~1)
- `info`: `:row` or `:col`

# Returns
`EpidemicBounds` with tighter lower bounds than the general case.

Based on Theorem 3.2 of Britton et al. (2025).
"""
function r0_bounds_detailed_balance(K::AbstractMatrix{<:Real}, π::AbstractVector{<:Real};
                                    info::Symbol=:row)
    info in (:row, :col, :both) || throw(ArgumentError("info must be :row, :col, or :both"))
    k = length(π)
    size(K) == (k, k) || throw(DimensionMismatch("K must be $k × $k"))
    all(x -> isfinite(x) && x >= 0, K) || throw(ArgumentError("K must be non-negative"))
    all(x -> isfinite(x) && x >= 0, π) || throw(ArgumentError("π must be non-negative"))

    # Filter out zero-population groups
    active = findall(x -> x > 0, π)
    isempty(active) && return EpidemicBounds(0.0, 0.0)

    row_sums = vec(sum(K; dims=2))
    col_sums = vec(sum(K; dims=1))
    π_active = π[active]

    if info == :row || info == :both
        r_active = row_sums[active]
        # r̄ = √(Σ πᵢ rᵢ²)
        r_bar = sqrt(sum(π_active[i] * r_active[i]^2 for i in eachindex(active)))
        r_upper = maximum(r_active)
        row_bounds = EpidemicBounds(r_bar, r_upper)
    end

    if info == :col || info == :both
        c_active = col_sums[active]
        π_inv_sum = sum(1.0 / π_active[i] for i in eachindex(active))
        # c̃ = √(Σ cⱼ²/πⱼ / Σ 1/πᵢ)
        c_tilde = sqrt(sum(c_active[i]^2 / π_active[i] for i in eachindex(active)) / π_inv_sum)
        c_upper = maximum(c_active)
        col_bounds = EpidemicBounds(c_tilde, c_upper)
    end

    if info == :row
        return row_bounds
    elseif info == :col
        return col_bounds
    else  # :both
        lo = max(row_bounds.lower, col_bounds.lower)
        hi = min(row_bounds.upper, col_bounds.upper)
        return EpidemicBounds(lo, hi)
    end
end

"""
    r0_bounds_detailed_balance(cm::ContactMatrix; info=:row, transmissibility=1, recovery_rate=1)

Convenience wrapper for `ContactMatrix` inputs.
"""
function r0_bounds_detailed_balance(cm::ContactMatrix; info::Symbol=:row,
                                    transmissibility::Real=1, recovery_rate::Real=1)
    K_contact = next_generation_matrix(cm;
        transmissibility=transmissibility, recovery_rate=recovery_rate)
    K_paper = Matrix(transpose(K_contact))
    pop = population(cm)
    π = pop ./ sum(pop)
    r0_bounds_detailed_balance(K_paper, π; info=info)
end

# ---------------------------------------------------------------------------
# Final size bounds
# ---------------------------------------------------------------------------

# κ(x) = max(0, x - 1 - log(x)), auxiliary for final-size bounds
function κ(x::Real)
    x = Float64(x)
    x <= 1.0 && return 0.0
    x - 1.0 - log(x)
end

"""
    final_size_bounds(K, π; info=:col)

Compute per-type final size bounds τᵢ given partial NGM information.

# Arguments
- `K`: Next-generation matrix (paper convention)
- `π`: Population fractions
- `info`: `:col` (Theorem 3.3) or `:row` (Theorem 3.4)

# Returns
`EpidemicBounds` where `lower` and `upper` are vectors of per-type bounds.
"""
function final_size_bounds(K::AbstractMatrix{<:Real}, π::AbstractVector{<:Real};
                           info::Symbol=:col)
    info in (:row, :col) || throw(ArgumentError("info must be :row or :col"))
    k = length(π)
    size(K) == (k, k) || throw(DimensionMismatch("K must be $k × $k"))

    active = findall(x -> x > 0, π)
    lower = zeros(k)
    upper = zeros(k)

    if info == :col
        # Theorem 3.3: known column sums
        col_sums = vec(sum(K; dims=1))

        # y* = min_j(πⱼ · t_{cⱼ}), y* = max_j(πⱼ · t_{cⱼ})
        yt_vals = [π[j] * solve_final_size_scalar(col_sums[j]) for j in active]
        y_lower = minimum(yt_vals)
        y_upper = maximum(yt_vals)

        for i in 1:k
            if π[i] == 0
                lower[i] = 0.0
                upper[i] = 0.0
            else
                lower[i] = 1.0 - exp(-col_sums[i] / π[i] * y_lower)
                upper[i] = 1.0 - exp(-col_sums[i] / π[i] * y_upper)
            end
        end
    else
        # Theorem 3.4: known row sums
        row_sums = vec(sum(K; dims=2))

        for i in 1:k
            if π[i] == 0
                lower[i] = 0.0
                upper[i] = 0.0
                continue
            end
            lower[i] = 0.0  # always valid

            # Kᵢ = (1/πᵢ) Σ_{j≠i} πⱼ κ(rⱼ)
            Ki = (1.0 / π[i]) * sum(π[j] * κ(row_sums[j]) for j in active if j != i;
                                    init=0.0)
            upper[i] = solve_final_size_ext(row_sums[i], Ki)
        end
    end

    EpidemicBounds(lower, upper)
end

"""
    final_size_bounds(cm::ContactMatrix; info=:col, transmissibility=1, recovery_rate=1)

Convenience wrapper for `ContactMatrix` inputs.
"""
function final_size_bounds(cm::ContactMatrix; info::Symbol=:col,
                           transmissibility::Real=1, recovery_rate::Real=1)
    K_contact = next_generation_matrix(cm;
        transmissibility=transmissibility, recovery_rate=recovery_rate)
    K_paper = Matrix(transpose(K_contact))
    pop = population(cm)
    π = pop ./ sum(pop)
    final_size_bounds(K_paper, π; info=info)
end

"""
    total_final_size_bounds(K, π; info=:row)

Compute bounds on the total final epidemic fraction τ̄ = Σᵢ πᵢτᵢ.

# Arguments
- `K`: Next-generation matrix (paper convention)
- `π`: Population fractions
- `info`: `:row` (Theorem 3.5) or `:col` (from Theorem 3.3)

# Returns
`EpidemicBounds` with scalar lower and upper bounds on τ̄.
"""
function total_final_size_bounds(K::AbstractMatrix{<:Real}, π::AbstractVector{<:Real};
                                 info::Symbol=:row)
    info in (:row, :col) || throw(ArgumentError("info must be :row or :col"))
    k = length(π)
    size(K) == (k, k) || throw(DimensionMismatch("K must be $k × $k"))

    active = findall(x -> x > 0, π)

    if info == :col
        # From Theorem 3.3: weighted sum of per-type bounds
        bounds = final_size_bounds(K, π; info=:col)
        lo = sum(π[i] * bounds.lower[i] for i in active; init=0.0)
        hi = sum(π[i] * bounds.upper[i] for i in active; init=0.0)
        return EpidemicBounds(lo, hi)
    end

    # Theorem 3.5: known row sums
    row_sums = vec(sum(K; dims=2))

    # Lower bound: min_i(πᵢ · t_{rᵢ})
    lo_vals = [π[i] * solve_final_size_scalar(row_sums[i]) for i in active]
    lo = minimum(lo_vals)

    # Upper bound τ*: solve for λ* via Eq. (11)
    # Check if all rᵢ ≤ 1 → τ* = 0
    if all(row_sums[i] <= 1.0 for i in active)
        return EpidemicBounds(lo, 0.0)
    end

    hi = _solve_total_final_size_upper(row_sums, π, active)
    return EpidemicBounds(lo, hi)
end

# Solve for the upper bound τ* on total final size (Theorem 3.5).
function _solve_total_final_size_upper(row_sums::Vector{Float64},
                                        π::AbstractVector{<:Real},
                                        active::Vector{Int})
    # The equation for λ*: LHS(λ) = RHS(λ)
    # LHS(λ) = Σᵢ πᵢ [log((1 + λrᵢ)/λ)]₊
    # RHS(λ) = Σᵢ rᵢπᵢ [1 - λ/(1 + λrᵢ)]₊
    # Equivalently: h(λ) = LHS(λ) - RHS(λ) = 0

    function h(λ)
        lhs = 0.0
        rhs = 0.0
        for i in active
            ri = row_sums[i]
            πi = Float64(π[i])
            val = log((1.0 + λ * ri) / λ)
            if val > 0
                lhs += πi * val
            end
            frac = 1.0 - λ / (1.0 + λ * ri)
            if frac > 0
                rhs += ri * πi * frac
            end
        end
        lhs - rhs
    end

    # Bisection on λ ∈ (0, large)
    # As λ → 0⁺: LHS → ∞, RHS → Σ rᵢπᵢ (bounded), so h → +∞
    # As λ → ∞: log((1+λr)/λ) = log(r + 1/λ) → log(r), bounded;
    #            1 - λ/(1+λr) = 1/(1+λr) → 0, so RHS → 0
    # Need to find where h crosses zero
    λ_lo = 1e-10
    λ_hi = 100.0

    # Expand upper bound if needed
    while h(λ_hi) > 0 && λ_hi < 1e10
        λ_hi *= 10.0
    end

    if h(λ_hi) > 0
        # h is always positive — use fallback
        return sum(Float64(π[i]) * solve_final_size_scalar(row_sums[i]) for i in active)
    end

    # Bisection
    for _ in 1:200
        λ_mid = (λ_lo + λ_hi) / 2
        if h(λ_mid) > 0
            λ_lo = λ_mid
        else
            λ_hi = λ_mid
        end
        (λ_hi - λ_lo) < 1e-12 * λ_lo && break
    end
    λ_star = (λ_lo + λ_hi) / 2

    # Compute τ* = Σᵢ πᵢ [1 - λ*/(1 + λ*rᵢ)]₊
    τ_star = 0.0
    for i in active
        ri = row_sums[i]
        val = 1.0 - λ_star / (1.0 + λ_star * ri)
        if val > 0
            τ_star += Float64(π[i]) * val
        end
    end
    τ_star
end

"""
    total_final_size_bounds(cm::ContactMatrix; info=:row, transmissibility=1, recovery_rate=1)

Convenience wrapper for `ContactMatrix` inputs.
"""
function total_final_size_bounds(cm::ContactMatrix; info::Symbol=:row,
                                 transmissibility::Real=1, recovery_rate::Real=1)
    K_contact = next_generation_matrix(cm;
        transmissibility=transmissibility, recovery_rate=recovery_rate)
    K_paper = Matrix(transpose(K_contact))
    pop = population(cm)
    π = pop ./ sum(pop)
    total_final_size_bounds(K_paper, π; info=info)
end

# ---------------------------------------------------------------------------
# Fiber uncertainty: R₀ and final-size ranges over the q-parameter family
# ---------------------------------------------------------------------------

"""
    epidemic_uncertainty(matrices; transmissibility=1, recovery_rate=1)

Compute the range of R₀ and total final size across a collection of contact
matrices (e.g., from sampling or MCMC over the q-parameter fiber).

# Returns
Named tuple with fields `r0::EpidemicBounds` and `final_size::EpidemicBounds`.
"""
function epidemic_uncertainty(matrices::AbstractVector{<:ContactMatrix};
                              transmissibility::Real=1, recovery_rate::Real=1)
    isempty(matrices) && throw(ArgumentError("matrices must be non-empty"))

    r0_vals = Float64[]
    fs_vals = Float64[]

    for cm in matrices
        K_contact = next_generation_matrix(cm;
            transmissibility=transmissibility, recovery_rate=recovery_rate)
        K_paper = Matrix(transpose(K_contact))
        pop = population(cm)
        π = pop ./ sum(pop)
        push!(r0_vals, maximum(abs.(eigvals(K_paper))))

        τ = solve_final_size_vector(K_paper, π)
        τ_bar = sum(π[i] * τ[i] for i in eachindex(π) if π[i] > 0; init=0.0)
        push!(fs_vals, τ_bar)
    end

    (r0=EpidemicBounds(minimum(r0_vals), maximum(r0_vals)),
     final_size=EpidemicBounds(minimum(fs_vals), maximum(fs_vals)))
end
