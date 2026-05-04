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

BlockAssortativityParams(; q::Dict{Symbol,Float64}=Dict{Symbol,Float64}()) =
    BlockAssortativityParams(q)

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
                block .= max.(block, 0.0)
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

# Apply q-parameter assortativity adjustment to a block.
# For each dimension with a q parameter, entries where the row and column groups
# share the same factor value in that dimension receive mass shifted from
# cross-group entries. The shift preserves row marginals exactly.
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

    block
end
