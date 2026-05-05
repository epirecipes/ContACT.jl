"""
Constrained generalized lifts for partial-data reconstruction.

This layer supports CoMix/SEP-style workflows where a rectangular
source-stratified intermediate matrix constrains reconstruction of a full square
matrix over the source product partition.
"""

"""
    ConstrainedGeneralizedLift(intermediate; source_map, structural_zeros=nothing)
    ConstrainedGeneralizedLift(full_partition, intermediate; source_map, structural_zeros=nothing)

Specification for reconstructing a full generalized contact matrix from partial
source-stratified observations.

`intermediate` has target rows over a base partition `P` and source columns over
a richer partition `S`, typically `P × Q`. `source_map` maps `S → P`, usually a
product projection. `structural_zeros` marks impossible source groups whose full
rows/columns must remain zero; if omitted, zero-population source groups are
treated as structural zeros.
"""
struct ConstrainedGeneralizedLift{P<:AbstractPartition,I<:SourceStratifiedContactMatrix}
    full_partition::P
    intermediate::I
    source_map::PartitionMap
    structural_zeros::Vector{Bool}
end

function ConstrainedGeneralizedLift(
    full_partition::P,
    intermediate::I;
    source_map::PartitionMap=PartitionMap(full_partition, target_partition(intermediate)),
    structural_zeros=nothing,
) where {P<:AbstractPartition,I<:SourceStratifiedContactMatrix}
    same_partition(full_partition, source_partition(intermediate)) || throw(ArgumentError(
        "full partition must match the intermediate source partition"))
    same_partition(source_map.domain, full_partition) || throw(ArgumentError(
        "source_map domain must match the full partition"))
    same_partition(source_map.codomain, target_partition(intermediate)) || throw(ArgumentError(
        "source_map codomain must match the intermediate target partition"))

    zeros_mask = if structural_zeros === nothing
        population(intermediate) .== 0
    else
        mask = Bool.(structural_zeros)
        length(mask) == n_groups(full_partition) || throw(DimensionMismatch(
            "structural_zeros length $(length(mask)) does not match full partition groups $(n_groups(full_partition))"))
        mask
    end

    counts = source_total_contacts(intermediate)
    for source in eachindex(zeros_mask)
        if zeros_mask[source] && !all(iszero, counts[:, source])
            throw(ArgumentError(
                "structural-zero source group $source has nonzero intermediate contacts"))
        end
    end

    ConstrainedGeneralizedLift{P,I}(full_partition, intermediate, source_map, zeros_mask)
end

function ConstrainedGeneralizedLift(intermediate::SourceStratifiedContactMatrix; kwargs...)
    ConstrainedGeneralizedLift(source_partition(intermediate), intermediate; kwargs...)
end

"""Extract the full partition for a constrained generalized lift."""
full_partition(spec::ConstrainedGeneralizedLift) = spec.full_partition

"""Extract the intermediate source-stratified matrix for a constrained lift."""
intermediate_matrix(spec::ConstrainedGeneralizedLift) = spec.intermediate

"""Extract the source-to-base projection used by a constrained lift."""
source_map(spec::ConstrainedGeneralizedLift) = spec.source_map

"""Return the structural-zero mask for a constrained lift."""
structural_zeros(spec::ConstrainedGeneralizedLift) = spec.structural_zeros

"""
    constrained_generalize(base, spec::ConstrainedGeneralizedLift)

Reconstruct a full generalized `ContactMatrix` from a reciprocal base matrix and
a source-stratified intermediate matrix.

The current solver uses a deterministic proportionate block transport. For each
base block, it preserves the aligned intermediate source marginals, enforces
reciprocal total contacts, and leaves structural-zero groups with zero rows and
columns. The base matrix must be reciprocal in total-contact space; the
alignment step checks this explicitly, which guarantees matching block totals for
off-diagonal transport plans.
"""
function constrained_generalize(base::ContactMatrix, spec::ConstrainedGeneralizedLift)
    aligned = align_source_stratified_matrix(
        spec.intermediate,
        base,
        spec.source_map,
    )

    full_partition = spec.full_partition
    full_pop = population(aligned)
    source_base = collect(spec.source_map.mapping)
    n_base = n_groups(base.partition)
    n_full = n_groups(full_partition)
    intermediate_counts = source_total_contacts(aligned)
    full_counts = zeros(Float64, n_full, n_full)

    for target_base in 1:n_base
        target_groups = findall(==(target_base), source_base)
        for source_base_group in target_base:n_base
            source_groups = findall(==(source_base_group), source_base)

            if target_base == source_base_group
                marginal = [intermediate_counts[target_base, g] for g in target_groups]
                block = _proportionate_transport(marginal, marginal)
                for (i, target_group) in enumerate(target_groups)
                    for (j, source_group) in enumerate(source_groups)
                        full_counts[target_group, source_group] = block[i, j]
                    end
                end
            else
                row_marginal = [intermediate_counts[source_base_group, g]
                                for g in target_groups]
                col_marginal = [intermediate_counts[target_base, g]
                                for g in source_groups]
                block = _proportionate_transport(row_marginal, col_marginal)
                for (i, target_group) in enumerate(target_groups)
                    for (j, source_group) in enumerate(source_groups)
                        value = block[i, j]
                        full_counts[target_group, source_group] = value
                        full_counts[source_group, target_group] = value
                    end
                end
            end
        end
    end

    for group in 1:n_full
        if spec.structural_zeros[group]
            if !all(iszero, full_counts[:, group]) || !all(iszero, full_counts[group, :])
                throw(ArgumentError("structural-zero group $group received reconstructed contacts"))
            end
        end
    end

    full_matrix = _source_mean_contacts_from_counts(full_counts, full_pop)
    ContactMatrix(full_matrix, full_partition, full_pop, MeanContacts())
end

generalize(base::ContactMatrix, spec::ConstrainedGeneralizedLift) =
    constrained_generalize(base, spec)

# ---------------------------------------------------------------------------
# Parameterized (q-parameter) block solver
# ---------------------------------------------------------------------------

"""
    BlockAssortativityParams(; q=Dict{Symbol,Float64}())

Per-block assortativity parameters. Each entry maps a product dimension symbol
(e.g. `:sep` or `:edu`) to a scalar `q` value in `[-1, 1]` controlling
mixing preference within that dimension:

- `q = 0`: proportionate mixing (default)
- `q > 0`: assortative mixing (same-group preference)
- `q < 0`: disassortative mixing (cross-group preference)

The parametrization adjusts the proportionate transport plan within each base
block by shifting mass toward or away from same-dimension-value entries while
preserving row and column marginals.
"""
struct BlockAssortativityParams
    q::Dict{Symbol,Float64}
end

function BlockAssortativityParams(; q::AbstractDict{Symbol,<:Real}=Dict{Symbol,Float64}())
    vals = Dict{Symbol,Float64}()
    for (dim, val) in q
        qval = Float64(val)
        isfinite(qval) || throw(ArgumentError("q-parameter for $dim must be finite"))
        abs(qval) <= 1.0 || throw(ArgumentError(
            "q-parameter for $dim must lie in [-1, 1], got $qval"))
        vals[dim] = qval
    end
    BlockAssortativityParams(vals)
end

"""
    ParameterizedConstrainedLift(spec; block_params=nothing, default_params=BlockAssortativityParams())

Wrap a `ConstrainedGeneralizedLift` with per-block or global assortativity
parameters for the q-parameter solver.

`block_params` is an optional dictionary mapping `(base_row, base_col)` tuples
to `BlockAssortativityParams`. If a block is not in `block_params`, the
`default_params` apply.
"""
struct ParameterizedConstrainedLift{S<:ConstrainedGeneralizedLift}
    spec::S
    block_params::Dict{Tuple{Int,Int},BlockAssortativityParams}
    default_params::BlockAssortativityParams
end

function ParameterizedConstrainedLift(spec::ConstrainedGeneralizedLift;
    block_params::Union{Nothing,Dict{Tuple{Int,Int},BlockAssortativityParams}}=nothing,
    default_params::BlockAssortativityParams=BlockAssortativityParams())
    bp = block_params === nothing ? Dict{Tuple{Int,Int},BlockAssortativityParams}() : block_params
    ParameterizedConstrainedLift(spec, bp, default_params)
end

"""
    constrained_generalize(base, pspec::ParameterizedConstrainedLift)

Reconstruct a full generalized `ContactMatrix` using per-block q-parameter
assortativity. Each base block is solved independently: the proportionate
transport plan is adjusted by shifting mass toward same-dimension-value entries
(assortative, q > 0) or away (disassortative, q < 0).

Throws `ArgumentError` if the resulting block has negative entries (infeasible
parameters).
"""
function constrained_generalize(base::ContactMatrix, pspec::ParameterizedConstrainedLift)
    spec = pspec.spec
    aligned = align_source_stratified_matrix(
        spec.intermediate,
        base,
        spec.source_map,
    )

    full_part = spec.full_partition
    full_pop = population(aligned)
    source_base = collect(spec.source_map.mapping)
    n_base = n_groups(base.partition)
    n_full = n_groups(full_part)
    intermediate_counts = source_total_contacts(aligned)
    full_counts = zeros(Float64, n_full, n_full)

    for target_base in 1:n_base
        target_groups = findall(==(target_base), source_base)
        for source_base_group in target_base:n_base
            source_groups = findall(==(source_base_group), source_base)
            block_key = (target_base, source_base_group)
            params = get(pspec.block_params, block_key, pspec.default_params)

            if target_base == source_base_group
                marginal = [intermediate_counts[target_base, g] for g in target_groups]
                block = _parameterized_transport(
                    marginal, marginal, target_groups, target_groups,
                    full_part, params)
                # Enforce symmetry for diagonal blocks (reciprocity)
                block = (block + transpose(block)) / 2
                block = _balance_transport(block, marginal, marginal)
                for (i, tg) in enumerate(target_groups)
                    for (j, sg) in enumerate(target_groups)
                        full_counts[tg, sg] = block[i, j]
                    end
                end
            else
                row_marginal = [intermediate_counts[source_base_group, g]
                                for g in target_groups]
                col_marginal = [intermediate_counts[target_base, g]
                                for g in source_groups]
                block = _parameterized_transport(
                    row_marginal, col_marginal, target_groups, source_groups,
                    full_part, params)
                for (i, tg) in enumerate(target_groups)
                    for (j, sg) in enumerate(source_groups)
                        value = block[i, j]
                        full_counts[tg, sg] = value
                        full_counts[sg, tg] = value
                    end
                end
            end
        end
    end

    # Final symmetrization of total contacts (Sinkhorn introduces ~1e-10 rounding)
    full_counts .= (full_counts .+ full_counts') ./ 2

    for group in 1:n_full
        if spec.structural_zeros[group]
            if !all(iszero, full_counts[:, group]) || !all(iszero, full_counts[group, :])
                throw(ArgumentError("structural-zero group $group received reconstructed contacts"))
            end
        end
    end

    full_matrix = _source_mean_contacts_from_counts(full_counts, full_pop)
    ContactMatrix(full_matrix, full_part, full_pop, MeanContacts())
end

generalize(base::ContactMatrix, pspec::ParameterizedConstrainedLift) =
    constrained_generalize(base, pspec)

"""
    is_feasible(base, pspec::ParameterizedConstrainedLift)

Check whether the parameterized constrained lift produces a non-negative matrix.
Returns `true` if all entries are non-negative, `false` otherwise.
"""
function is_feasible(base::ContactMatrix, pspec::ParameterizedConstrainedLift)
    try
        cm = constrained_generalize(base, pspec)
        return all(matrix(cm) .>= -1e-12)
    catch e
        e isa ArgumentError && return false
        rethrow(e)
    end
end

"""
    sample_constrained_lifts(base, spec::ConstrainedGeneralizedLift, n::Int;
        dimensions=Symbol[], bounds=(-1.0, 1.0), rng=Random.default_rng())

Sample `n` feasible parameterized lifts by uniformly sampling q-parameters
within `bounds` and rejecting infeasible (negative-entry) solutions.

Returns a vector of `(params, matrix)` tuples for accepted samples.
"""
function sample_constrained_lifts(base::ContactMatrix, spec::ConstrainedGeneralizedLift, n::Int;
    dimensions::AbstractVector{Symbol}=Symbol[],
    bounds::Tuple{Real,Real}=(-1.0, 1.0),
    rng=nothing)
    n > 0 || throw(ArgumentError("n must be positive"))
    bounds[1] < bounds[2] || throw(ArgumentError("bounds must satisfy lower < upper"))

    dims = if isempty(dimensions)
        full_part = spec.full_partition
        full_part isa ProductPartition || throw(ArgumentError(
            "cannot infer dimensions from non-product partition; specify dimensions explicitly"))
        base_dim = dimension(spec.source_map.codomain)
        [dimension(f) for f in full_part.factors if dimension(f) != base_dim]
    else
        collect(dimensions)
    end

    _rng = rng === nothing ? Random.default_rng() : rng
    lo, hi = Float64(bounds[1]), Float64(bounds[2])
    results = Vector{Tuple{BlockAssortativityParams,ContactMatrix}}()

    attempts = 0
    max_attempts = n * 100

    while length(results) < n && attempts < max_attempts
        attempts += 1
        q_vals = Dict(d => lo + (hi - lo) * rand(_rng) for d in dims)
        params = BlockAssortativityParams(q=q_vals)
        pspec = ParameterizedConstrainedLift(spec; default_params=params)
        if is_feasible(base, pspec)
            cm = constrained_generalize(base, pspec)
            push!(results, (params, cm))
        end
    end

    results
end

# ---------------------------------------------------------------------------
# Per-block per-dimension q-parameter space and MCMC sampling
# ---------------------------------------------------------------------------

"""
    QParameterSpace(base, spec; dimensions=Symbol[], bounds=(-1.0, 1.0))

Defines the full per-block × per-dimension q-parameter space for a constrained
lift. Each non-trivial base block `(i,j)` with `i ≤ j` gets independent
q-parameters for each specified dimension.

The space is represented as a flat vector of parameters with named indexing.
"""
struct QParameterSpace
    block_keys::Vector{Tuple{Int,Int}}
    dimensions::Vector{Symbol}
    bounds::Tuple{Float64,Float64}
    n_params::Int
end

function QParameterSpace(base::ContactMatrix, spec::ConstrainedGeneralizedLift;
    dimensions::AbstractVector{Symbol}=Symbol[],
    bounds::Tuple{Real,Real}=(-1.0, 1.0))
    bounds[1] < bounds[2] || throw(ArgumentError("bounds must satisfy lower < upper"))

    dims = if isempty(dimensions)
        full_part = spec.full_partition
        full_part isa ProductPartition || throw(ArgumentError(
            "cannot infer dimensions from non-product partition; specify dimensions explicitly"))
        base_dim = dimension(spec.source_map.codomain)
        Symbol[dimension(f) for f in full_part.factors if dimension(f) != base_dim]
    else
        Symbol.(collect(dimensions))
    end

    n_base = n_groups(base.partition)
    block_keys = Tuple{Int,Int}[]
    for i in 1:n_base
        for j in i:n_base
            push!(block_keys, (i, j))
        end
    end

    QParameterSpace(block_keys, dims, (Float64(bounds[1]), Float64(bounds[2])),
                    length(block_keys) * length(dims))
end

function _vector_to_block_params(space::QParameterSpace, θ::AbstractVector{Float64})
    length(θ) == space.n_params || throw(DimensionMismatch(
        "parameter vector length $(length(θ)) ≠ expected $(space.n_params)"))
    block_params = Dict{Tuple{Int,Int},BlockAssortativityParams}()
    idx = 1
    for block_key in space.block_keys
        q_vals = Dict{Symbol,Float64}()
        for dim in space.dimensions
            q_vals[dim] = θ[idx]
            idx += 1
        end
        block_params[block_key] = BlockAssortativityParams(q=q_vals)
    end
    block_params
end

function _block_params_to_vector(space::QParameterSpace,
    block_params::Dict{Tuple{Int,Int},BlockAssortativityParams})
    θ = zeros(Float64, space.n_params)
    idx = 1
    for block_key in space.block_keys
        params = get(block_params, block_key, BlockAssortativityParams())
        for dim in space.dimensions
            θ[idx] = get(params.q, dim, 0.0)
            idx += 1
        end
    end
    θ
end

"""
    sample_perblock_lifts(base, spec, n; dimensions=Symbol[], bounds=(-1.0, 1.0), rng=nothing)

Sample `n` feasible parameterized lifts with independent per-block per-dimension
q-parameters. Uses rejection sampling over the full parameter space.

Returns a vector of `(block_params, matrix)` tuples.
"""
function sample_perblock_lifts(base::ContactMatrix, spec::ConstrainedGeneralizedLift, n::Int;
    dimensions::AbstractVector{Symbol}=Symbol[],
    bounds::Tuple{Real,Real}=(-1.0, 1.0),
    rng=nothing)
    n > 0 || throw(ArgumentError("n must be positive"))
    space = QParameterSpace(base, spec; dimensions=dimensions, bounds=bounds)
    _rng = rng === nothing ? Random.default_rng() : rng
    lo, hi = space.bounds

    results = Vector{Tuple{Dict{Tuple{Int,Int},BlockAssortativityParams},ContactMatrix}}()
    attempts = 0
    max_attempts = n * 200

    while length(results) < n && attempts < max_attempts
        attempts += 1
        θ = [lo + (hi - lo) * rand(_rng) for _ in 1:space.n_params]
        block_params = _vector_to_block_params(space, θ)
        pspec = ParameterizedConstrainedLift(spec; block_params=block_params)
        if is_feasible(base, pspec)
            cm = constrained_generalize(base, pspec)
            push!(results, (block_params, cm))
        end
    end

    results
end

"""
    MCMCResult

Result of MCMC sampling over the q-parameter space.

Fields:
- `chain`: Vector of accepted parameter vectors (each a `Dict{Tuple{Int,Int},BlockAssortativityParams}`)
- `matrices`: Corresponding contact matrices
- `log_densities`: Log-density at each sample
- `acceptance_rate`: Fraction of proposals accepted
- `space`: The `QParameterSpace` defining the parameter layout
"""
struct MCMCResult
    chain::Vector{Dict{Tuple{Int,Int},BlockAssortativityParams}}
    matrices::Vector{ContactMatrix}
    log_densities::Vector{Float64}
    acceptance_rate::Float64
    space::QParameterSpace
end

"""
    mcmc_constrained_lifts(base, spec, n_samples;
        dimensions=Symbol[], bounds=(-1.0, 1.0),
        log_density=nothing, proposal_scale=0.1,
        burnin=100, thin=1, rng=nothing)

Sample from the space of feasible q-parameterized contact matrices using
Metropolis-Hastings MCMC.

# Arguments
- `base`: Reciprocal base `ContactMatrix` over the coarse partition.
- `spec`: `ConstrainedGeneralizedLift` specification.
- `n_samples`: Number of post-burnin samples to collect.
- `dimensions`: Which product dimensions get q-parameters (auto-detected if empty).
- `bounds`: Hard bounds on q-parameters (proposals outside are rejected).
- `log_density`: Optional function `(ContactMatrix, block_params) -> Float64` giving
  the log-density to target. If `nothing`, samples uniformly over feasible region.
- `proposal_scale`: Standard deviation of Gaussian random walk proposals.
- `burnin`: Number of initial samples to discard.
- `thin`: Keep every `thin`-th sample after burnin.
- `rng`: Random number generator.

# Returns
An `MCMCResult` with the chain, matrices, and diagnostics.
"""
function mcmc_constrained_lifts(base::ContactMatrix, spec::ConstrainedGeneralizedLift,
    n_samples::Int;
    dimensions::AbstractVector{Symbol}=Symbol[],
    bounds::Tuple{Real,Real}=(-1.0, 1.0),
    log_density=nothing,
    proposal_scale::Real=0.1,
    burnin::Int=100,
    thin::Int=1,
    rng=nothing)

    n_samples > 0 || throw(ArgumentError("n_samples must be positive"))
    burnin >= 0 || throw(ArgumentError("burnin must be non-negative"))
    thin >= 1 || throw(ArgumentError("thin must be >= 1"))
    proposal_scale > 0 || throw(ArgumentError("proposal_scale must be positive"))

    space = QParameterSpace(base, spec; dimensions=dimensions, bounds=bounds)
    _rng = rng === nothing ? Random.default_rng() : rng
    lo, hi = space.bounds
    σ = Float64(proposal_scale)

    # Find a feasible starting point
    θ_current = _find_feasible_start(base, spec, space, _rng; max_attempts=1000)
    θ_current === nothing && throw(ArgumentError(
        "could not find a feasible starting point within bounds"))

    block_params_current = _vector_to_block_params(space, θ_current)
    pspec_current = ParameterizedConstrainedLift(spec; block_params=block_params_current)
    cm_current = constrained_generalize(base, pspec_current)
    ld_current = log_density === nothing ? 0.0 : Float64(log_density(cm_current, block_params_current))

    # MCMC loop
    total_steps = burnin + n_samples * thin
    chain = Vector{Dict{Tuple{Int,Int},BlockAssortativityParams}}()
    matrices = Vector{ContactMatrix}()
    log_densities = Vector{Float64}()
    accepted = 0

    for step in 1:total_steps
        # Propose: Gaussian random walk
        θ_proposal = θ_current .+ σ .* randn(_rng, space.n_params)

        # Check bounds
        if all(lo .<= θ_proposal .<= hi)
            block_params_prop = _vector_to_block_params(space, θ_proposal)
            pspec_prop = ParameterizedConstrainedLift(spec; block_params=block_params_prop)

            if is_feasible(base, pspec_prop)
                cm_prop = constrained_generalize(base, pspec_prop)
                ld_prop = log_density === nothing ? 0.0 :
                    Float64(log_density(cm_prop, block_params_prop))

                # Metropolis acceptance (symmetric proposal → ratio is just density ratio)
                log_α = ld_prop - ld_current
                if log_α >= 0 || log(rand(_rng)) < log_α
                    θ_current = θ_proposal
                    block_params_current = block_params_prop
                    cm_current = cm_prop
                    ld_current = ld_prop
                    accepted += 1
                end
            end
        end

        # Collect after burnin, respecting thin
        if step > burnin && (step - burnin) % thin == 0
            push!(chain, block_params_current)
            push!(matrices, cm_current)
            push!(log_densities, ld_current)
        end
    end

    MCMCResult(chain, matrices, log_densities, accepted / total_steps, space)
end

function _find_feasible_start(base::ContactMatrix, spec::ConstrainedGeneralizedLift,
    space::QParameterSpace, rng; max_attempts::Int=1000)
    lo, hi = space.bounds
    # Try q=0 first (always feasible if proportionate solver works)
    θ_zero = zeros(Float64, space.n_params)
    block_params = _vector_to_block_params(space, θ_zero)
    pspec = ParameterizedConstrainedLift(spec; block_params=block_params)
    is_feasible(base, pspec) && return θ_zero

    # Random search
    for _ in 1:max_attempts
        θ = [lo + (hi - lo) * rand(rng) for _ in 1:space.n_params]
        block_params = _vector_to_block_params(space, θ)
        pspec = ParameterizedConstrainedLift(spec; block_params=block_params)
        is_feasible(base, pspec) && return θ
    end
    nothing
end

function _proportionate_transport(row_marginal::AbstractVector{<:Real},
                                  col_marginal::AbstractVector{<:Real})
    row = Float64.(row_marginal)
    col = Float64.(col_marginal)
    all(x -> isfinite(x) && x >= 0, row) ||
        throw(ArgumentError("row marginals must be finite and non-negative"))
    all(x -> isfinite(x) && x >= 0, col) ||
        throw(ArgumentError("column marginals must be finite and non-negative"))

    row_total = sum(row)
    col_total = sum(col)
    isapprox(row_total, col_total; rtol=1e-8, atol=1e-8) || throw(ArgumentError(
        "block marginals must have equal totals, got $row_total and $col_total"))
    row_total == 0 && return zeros(Float64, length(row), length(col))
    row * transpose(col) ./ row_total
end

# Apply q-parameter assortativity adjustment to a block, then rebalance the
# candidate plan so both row and column marginals are preserved.
function _parameterized_transport(row_marginal::AbstractVector{<:Real},
                                  col_marginal::AbstractVector{<:Real},
                                  row_groups::Vector{Int},
                                  col_groups::Vector{Int},
                                  full_partition::AbstractPartition,
                                  params::BlockAssortativityParams)
    row = Float64.(row_marginal)
    col = Float64.(col_marginal)
    all(x -> isfinite(x) && x >= 0, row) ||
        throw(ArgumentError("row marginals must be finite and non-negative"))
    all(x -> isfinite(x) && x >= 0, col) ||
        throw(ArgumentError("column marginals must be finite and non-negative"))

    row_total = sum(row)
    col_total = sum(col)
    isapprox(row_total, col_total; rtol=1e-8, atol=1e-8) || throw(ArgumentError(
        "block marginals must have equal totals, got $row_total and $col_total"))
    row_total == 0 && return zeros(Float64, length(row), length(col))

    # Start with proportionate transport
    block = row * transpose(col) ./ row_total

    # Apply q-parameter adjustment for each dimension
    isempty(params.q) && return block
    full_partition isa ProductPartition || return block

    sizes = Tuple(n_groups(f) for f in full_partition.factors)

    for (dim, q) in params.q
        abs(q) < 1e-15 && continue
        factor_idx = findfirst(f -> dimension(f) == dim, full_partition.factors)
        factor_idx === nothing && continue

        nrow = length(row_groups)
        ncol = length(col_groups)

        for i in 1:nrow
            row_factor_val = _cartesian_indices(row_groups[i], sizes)[factor_idx]
            same_mass = 0.0
            diff_mass = 0.0
            for j in 1:ncol
                col_factor_val = _cartesian_indices(col_groups[j], sizes)[factor_idx]
                if row_factor_val == col_factor_val
                    same_mass += block[i, j]
                else
                    diff_mass += block[i, j]
                end
            end

            shift_budget = min(same_mass, diff_mass)
            (shift_budget > 0 && same_mass > 0 && diff_mass > 0) || continue
            shift_amount = q * shift_budget

            for j in 1:ncol
                col_factor_val = _cartesian_indices(col_groups[j], sizes)[factor_idx]
                if row_factor_val == col_factor_val
                    block[i, j] += shift_amount * (block[i, j] / same_mass)
                else
                    block[i, j] -= shift_amount * (block[i, j] / diff_mass)
                end
            end
        end

        if any(block .< -1e-10)
            throw(ArgumentError(
                "q-parameter q_$(dim) = $q produces negative entries; infeasible parameters"))
        end
        block .= max.(block, 0.0)
    end

    _balance_transport(block, row, col)
end

function _balance_transport(candidate::AbstractMatrix{<:Real},
                            row_marginal::AbstractVector{<:Real},
                            col_marginal::AbstractVector{<:Real};
                            rtol::Real=1e-10,
                            atol::Real=1e-10,
                            maxiter::Int=10_000)
    row = Float64.(row_marginal)
    col = Float64.(col_marginal)
    row_total = sum(row)
    col_total = sum(col)
    isapprox(row_total, col_total; rtol=1e-8, atol=1e-8) || throw(ArgumentError(
        "block marginals must have equal totals, got $row_total and $col_total"))
    row_total == 0 && return zeros(Float64, length(row), length(col))

    block = Float64.(candidate)
    size(block) == (length(row), length(col)) || throw(DimensionMismatch(
        "candidate block size $(size(block)) does not match marginals $(length(row)) × $(length(col))"))
    all(isfinite, block) || throw(ArgumentError("candidate transport entries must be finite"))
    any(block .< -atol) && throw(ArgumentError("candidate transport has negative entries"))
    block .= max.(block, 0.0)

    for i in eachindex(row)
        if row[i] > 0 && sum(block[i, :]) <= 0
            throw(ArgumentError("cannot balance transport: row $i has positive marginal but zero support"))
        end
    end
    for j in eachindex(col)
        if col[j] > 0 && sum(block[:, j]) <= 0
            throw(ArgumentError("cannot balance transport: column $j has positive marginal but zero support"))
        end
    end

    scale = max(row_total, 1.0)
    for _ in 1:maxiter
        for i in eachindex(row)
            current = sum(block[i, :])
            if row[i] == 0
                block[i, :] .= 0.0
            else
                current > 0 || throw(ArgumentError(
                    "cannot balance transport: row $i has positive marginal but zero support"))
                block[i, :] .*= row[i] / current
            end
        end
        for j in eachindex(col)
            current = sum(block[:, j])
            if col[j] == 0
                block[:, j] .= 0.0
            else
                current > 0 || throw(ArgumentError(
                    "cannot balance transport: column $j has positive marginal but zero support"))
                block[:, j] .*= col[j] / current
            end
        end

        row_err = maximum(abs.(vec(sum(block; dims=2)) .- row))
        col_err = maximum(abs.(vec(sum(block; dims=1)) .- col))
        max(row_err, col_err) <= max(atol, rtol * scale) && return block
    end

    throw(ArgumentError("could not balance q-parameter transport to requested marginals"))
end
