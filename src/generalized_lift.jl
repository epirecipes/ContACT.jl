"""
Generalized contact-matrix lifts.

These operations construct Manna-style generalized contact matrices by lifting a
matrix over a base partition `P` to a product partition `P × Q`. The lift is a
section of the canonical projection `P × Q → P`: coarsening the result back to
`P` recovers the input matrix exactly.
"""

"""
    GeneralizedMixingKernel

Abstract supertype for blockwise mixing assumptions used by `generalize`.

Rows of a block kernel index the target group in the added partition, and
columns index the source group. Each block kernel must be non-negative and sum
to one, so it splits one base total-contact block across `Q × Q`.
"""
abstract type GeneralizedMixingKernel end

"""
    RandomMixing([distribution])

Product/random mixing over the added partition. If `distribution` is omitted,
the marginal distribution is inferred from the product population supplied to
`GeneralizedLift`.
"""
struct RandomMixing <: GeneralizedMixingKernel
    distribution::Union{Nothing,Vector{Float64}}
end

RandomMixing() = RandomMixing(nothing)
RandomMixing(distribution::AbstractVector{<:Real}) =
    RandomMixing(_validated_simplex(distribution, "random mixing distribution"))

"""
    BlockMixing(diagonal, offdiagonal=diagonal)

Explicit block kernels for a generalized lift.

`diagonal` is used for base blocks `(i, i)` and must be symmetric to preserve
reciprocity within a base group. `offdiagonal` is used for lower-triangular base
blocks, while its transpose is used for the reverse direction.
"""
struct BlockMixing <: GeneralizedMixingKernel
    diagonal::Matrix{Float64}
    offdiagonal::Matrix{Float64}

    function BlockMixing(diagonal::AbstractMatrix{<:Real},
                         offdiagonal::AbstractMatrix{<:Real})
        diag = _validate_standalone_block(diagonal, "diagonal block"; symmetric=true)
        offdiag = _validate_standalone_block(offdiagonal, "offdiagonal block")
        size(diag) == size(offdiag) || throw(DimensionMismatch(
            "diagonal and offdiagonal blocks must have the same size"))
        new(diag, offdiag)
    end
end

BlockMixing(block::AbstractMatrix{<:Real}) = BlockMixing(block, block)

"""
    AssortativeDimensionMixing(activity, assortativity; offdiag_split=nothing)

Parameterized Manna-style mixing over a three-level added dimension.

`activity` is the share of contacts assigned to each added group and must sum to
one. `assortativity[k]` is the fraction of group `k`'s activity kept within
group `k`. For three groups, `offdiag_split = [d₁, d₂, d₃]` reproduces the
off-diagonal construction used by Manna et al.; if omitted, the same symmetric
assortative block is used for all base blocks.
"""
struct AssortativeDimensionMixing <: GeneralizedMixingKernel
    activity::Vector{Float64}
    assortativity::Vector{Float64}
    offdiag_split::Union{Nothing,Vector{Float64}}
end

function AssortativeDimensionMixing(activity::AbstractVector{<:Real},
                                    assortativity::AbstractVector{<:Real};
                                    offdiag_split=nothing)
    activity_vec = _validated_simplex(activity, "activity")
    r = Float64.(assortativity)
    length(r) == length(activity_vec) || throw(DimensionMismatch(
        "assortativity length $(length(r)) does not match activity length $(length(activity_vec))"))
    all(x -> isfinite(x) && 0 <= x <= 1, r) ||
        throw(ArgumentError("assortativity values must be finite and in [0, 1]"))

    split = if offdiag_split === nothing
        nothing
    else
        s = Float64.(offdiag_split)
        length(s) == length(activity_vec) || throw(DimensionMismatch(
            "offdiag_split length $(length(s)) does not match activity length $(length(activity_vec))"))
        all(x -> isfinite(x) && 0 <= x <= 1, s) ||
            throw(ArgumentError("offdiag_split values must be finite and in [0, 1]"))
        s
    end

    AssortativeDimensionMixing(activity_vec, r, split)
end

"""
    GeneralizedLift(partition; distribution, mixing=RandomMixing())
    GeneralizedLift(partition, population; mixing=RandomMixing())

Specification for lifting a contact matrix over `P` to `P × partition`.

Use `distribution` when the added dimension is independent of the base
partition. Use `population` when a full joint population over `P × partition` is
available; matrix-shaped populations are interpreted as `(base, added)`.
"""
struct GeneralizedLift{Q<:AbstractPartition,K<:GeneralizedMixingKernel}
    partition::Q
    population::Vector{Float64}
    population_is_distribution::Bool
    mixing::K
end

function GeneralizedLift(partition::Q, population::AbstractVector{<:Real};
                         mixing::Union{Symbol,GeneralizedMixingKernel}=RandomMixing()) where {Q<:AbstractPartition}
    pop = _validated_nonnegative_vector(population, "product population")
    GeneralizedLift{Q,typeof(_generalized_mixing(mixing))}(
        partition, pop, false, _generalized_mixing(mixing))
end

function GeneralizedLift(partition::Q, population::AbstractMatrix{<:Real};
                         mixing::Union{Symbol,GeneralizedMixingKernel}=RandomMixing()) where {Q<:AbstractPartition}
    size(population, 2) == n_groups(partition) || throw(DimensionMismatch(
        "population matrix must have one column per added partition group"))
    flat = Float64[population[i, q] for i in axes(population, 1) for q in axes(population, 2)]
    GeneralizedLift(partition, flat; mixing=mixing)
end

function GeneralizedLift(partition::Q;
                         distribution=nothing,
                         mixing::Union{Symbol,GeneralizedMixingKernel}=RandomMixing()) where {Q<:AbstractPartition}
    distribution === nothing &&
        throw(ArgumentError("GeneralizedLift(partition; ...) requires a distribution keyword"))
    dist = _validated_simplex(distribution, "added-dimension distribution")
    length(dist) == n_groups(partition) || throw(DimensionMismatch(
        "distribution length $(length(dist)) does not match added partition groups $(n_groups(partition))"))
    GeneralizedLift{Q,typeof(_generalized_mixing(mixing))}(
        partition, dist, true, _generalized_mixing(mixing))
end

"""
    product_population(base_population, added_distribution)

Create a row-major product-population vector for `P × Q` from a base population
and an independent distribution over `Q`.
"""
function product_population(base_population::AbstractVector{<:Real},
                            added_distribution::AbstractVector{<:Real})
    base = _validated_nonnegative_vector(base_population, "base population")
    dist = _validated_simplex(added_distribution, "added-dimension distribution")
    Float64[base[i] * dist[q] for i in eachindex(base) for q in eachindex(dist)]
end

"""
    generalize(cm, spec::GeneralizedLift)
    generalize(cm, partition; distribution, mixing=RandomMixing())
    generalize(cm, partition, population; mixing=RandomMixing())

Lift `cm` from partition `P` to `P × partition` using a blockwise mixing
assumption. The result coarsens exactly back to `cm` along the product
projection.
"""
function generalize(cm::ContactMatrix, spec::GeneralizedLift)
    cm.semantics isa MeanContacts || throw(ArgumentError(
        "generalized contact-matrix lifts require MeanContacts semantics"))

    base_partition = cm.partition
    added_partition = spec.partition
    n_base = n_groups(base_partition)
    n_added = n_groups(added_partition)
    sizes = (n_base, n_added)

    refined_partition = base_partition × added_partition
    refined_pop = _resolved_lift_population(cm, spec, n_base, n_added)
    diagonal_block, offdiagonal_block =
        _mixing_blocks(spec.mixing, n_base, n_added, refined_pop)

    base_counts = _total_contacts(cm)
    refined_counts = zeros(Float64, n_base * n_added, n_base * n_added)

    for target_base in 1:n_base
        for source_base in 1:n_base
            block = _lift_block(
                diagonal_block,
                offdiagonal_block,
                target_base,
                source_base,
            )
            total = base_counts[target_base, source_base]
            total == 0 && continue

            for target_added in 1:n_added
                target_idx = _cartesian_index((target_base, target_added), sizes)
                for source_added in 1:n_added
                    source_idx = _cartesian_index((source_base, source_added), sizes)
                    refined_counts[target_idx, source_idx] =
                        total * block[target_added, source_added]
                end
            end
        end
    end

    refined_matrix = _mean_contacts_from_total_counts(refined_counts, refined_pop)
    ContactMatrix(refined_matrix, refined_partition, refined_pop, MeanContacts())
end

generalize(cm::ContactMatrix, partition::AbstractPartition; kwargs...) =
    generalize(cm, GeneralizedLift(partition; kwargs...))

generalize(cm::ContactMatrix, partition::AbstractPartition,
           population::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}; kwargs...) =
    generalize(cm, GeneralizedLift(partition, population; kwargs...))

"""
    generalized_lift(args...; kwargs...)

Alias for `generalize`.
"""
generalized_lift(args...; kwargs...) = generalize(args...; kwargs...)

function _generalized_mixing(mixing::GeneralizedMixingKernel)
    mixing
end

function _generalized_mixing(mixing::Symbol)
    if mixing in (:random, :product, :proportionate, :proportional)
        RandomMixing()
    else
        throw(ArgumentError("unknown generalized mixing kernel :$mixing"))
    end
end

function _validated_nonnegative_vector(values::AbstractVector{<:Real}, name::AbstractString)
    v = Float64.(values)
    all(x -> isfinite(x) && x >= 0, v) ||
        throw(ArgumentError("$name values must be finite and non-negative"))
    v
end

function _validated_simplex(values::AbstractVector{<:Real}, name::AbstractString)
    v = _validated_nonnegative_vector(values, name)
    total = sum(v)
    total > 0 || throw(ArgumentError("$name must have positive total"))
    isapprox(total, 1.0; rtol=1e-8, atol=1e-10) ||
        throw(ArgumentError("$name must sum to 1, got $total"))
    v
end

function _resolved_lift_population(cm::ContactMatrix, spec::GeneralizedLift,
                                   n_base::Int, n_added::Int)
    refined_pop = if spec.population_is_distribution
        product_population(population(cm), spec.population)
    else
        spec.population
    end

    length(refined_pop) == n_base * n_added || throw(DimensionMismatch(
        "product population length $(length(refined_pop)) does not match $n_base × $n_added groups"))

    base_marginal = zeros(Float64, n_base)
    for base in 1:n_base
        for added in 1:n_added
            base_marginal[base] += refined_pop[_cartesian_index((base, added), (n_base, n_added))]
        end
    end
    isapprox(base_marginal, population(cm); rtol=1e-8, atol=1e-8) ||
        throw(ArgumentError("product population must marginalize to the ContactMatrix population"))

    refined_pop
end

function _total_contacts(cm::ContactMatrix)
    M = matrix(cm)
    pop = population(cm)
    counts = zeros(Float64, n_groups(cm), n_groups(cm))
    for source in 1:n_groups(cm)
        counts[:, source] .= M[:, source] .* pop[source]
    end
    counts
end

function _mean_contacts_from_total_counts(counts::AbstractMatrix{<:Real},
                                          pop::AbstractVector{<:Real})
    M = zeros(Float64, size(counts))
    for source in eachindex(pop)
        if pop[source] == 0
            all(iszero, counts[:, source]) || throw(ArgumentError(
                "cannot convert total contacts for zero-population product group $source"))
        else
            M[:, source] .= counts[:, source] ./ pop[source]
        end
    end
    M
end

function _mixing_blocks(mixing::RandomMixing, n_base::Int, n_added::Int,
                        refined_pop::AbstractVector{<:Real})
    dist = if mixing.distribution === nothing
        marginal = zeros(Float64, n_added)
        total = sum(refined_pop)
        total > 0 || throw(ArgumentError("cannot infer random mixing distribution from zero population"))
        for base in 1:n_base
            for added in 1:n_added
                marginal[added] += refined_pop[_cartesian_index((base, added), (n_base, n_added))]
            end
        end
        marginal ./ total
    else
        length(mixing.distribution) == n_added || throw(DimensionMismatch(
            "random mixing distribution length $(length(mixing.distribution)) does not match added partition groups $n_added"))
        mixing.distribution
    end

    block = dist * transpose(dist)
    _validate_block_matrix(block, n_added, "random mixing block"; symmetric=true)
    block, block
end

function _mixing_blocks(mixing::BlockMixing, ::Int, n_added::Int,
                        ::AbstractVector{<:Real})
    diagonal = _validate_block_matrix(mixing.diagonal, n_added, "diagonal block"; symmetric=true)
    offdiagonal = _validate_block_matrix(mixing.offdiagonal, n_added, "offdiagonal block")
    diagonal, offdiagonal
end

function _mixing_blocks(mixing::AssortativeDimensionMixing, ::Int, n_added::Int,
                        ::AbstractVector{<:Real})
    length(mixing.activity) == n_added || throw(DimensionMismatch(
        "activity length $(length(mixing.activity)) does not match added partition groups $n_added"))
    diagonal = _assortative_diagonal_block(mixing.activity, mixing.assortativity)
    offdiagonal = if mixing.offdiag_split === nothing
        diagonal
    else
        _assortative_offdiagonal_block(
            mixing.activity,
            mixing.assortativity,
            mixing.offdiag_split,
        )
    end

    diagonal = _validate_block_matrix(diagonal, n_added, "assortative diagonal block"; symmetric=true)
    offdiagonal = _validate_block_matrix(offdiagonal, n_added, "assortative offdiagonal block")
    diagonal, offdiagonal
end

function _validate_block_matrix(block::AbstractMatrix{<:Real}, n::Int, name::AbstractString;
                                symmetric::Bool=false)
    size(block) == (n, n) || throw(DimensionMismatch(
        "$name size $(size(block)) does not match $n added groups"))
    _validate_standalone_block(block, name; symmetric=symmetric)
end

function _validate_standalone_block(block::AbstractMatrix{<:Real}, name::AbstractString;
                                    symmetric::Bool=false)
    size(block, 1) == size(block, 2) || throw(DimensionMismatch(
        "$name must be square, got size $(size(block))"))
    M = Matrix{Float64}(block)
    all(x -> isfinite(x) && x >= 0, M) ||
        throw(ArgumentError("$name entries must be finite and non-negative"))
    total = sum(M)
    isapprox(total, 1.0; rtol=1e-8, atol=1e-10) ||
        throw(ArgumentError("$name must sum to 1, got $total"))
    if symmetric && !isapprox(M, transpose(M); rtol=1e-8, atol=1e-10)
        throw(ArgumentError("$name must be symmetric"))
    end
    M
end

function _assortative_diagonal_block(activity::Vector{Float64}, r::Vector{Float64})
    n = length(activity)
    M = zeros(Float64, n, n)
    if n == 1
        M[1, 1] = 1.0
        return M
    elseif n == 2
        residual = (1 .- r) .* activity
        isapprox(residual[1], residual[2]; rtol=1e-8, atol=1e-10) ||
            throw(ArgumentError(
                "two-group assortative diagonal block requires equal residual activity; use BlockMixing for a custom kernel"))
        M[1, 1] = r[1] * activity[1]
        M[2, 2] = r[2] * activity[2]
        M[1, 2] = residual[1]
        M[2, 1] = residual[2]
        return M
    elseif n == 3
        residual = (1 .- r) .* activity
        M[1, 1] = r[1] * activity[1]
        M[2, 2] = r[2] * activity[2]
        M[3, 3] = r[3] * activity[3]
        M[1, 2] = M[2, 1] = (residual[1] + residual[2] - residual[3]) / 2
        M[1, 3] = M[3, 1] = (residual[1] + residual[3] - residual[2]) / 2
        M[2, 3] = M[3, 2] = (residual[2] + residual[3] - residual[1]) / 2
        if any(<(-1e-10), M)
            throw(ArgumentError(
                "assortative diagonal parameters produce negative off-diagonal mass; use BlockMixing for a custom kernel"))
        end
        return max.(M, 0.0)
    end

    throw(ArgumentError(
        "automatic assortative diagonal blocks are implemented for one, two, or three groups; use BlockMixing for $n groups"))
end

function _assortative_offdiagonal_block(activity::Vector{Float64}, r::Vector{Float64},
                                        split::Vector{Float64})
    n = length(activity)
    n == 3 || throw(ArgumentError(
        "offdiag_split currently follows the three-group Manna et al. construction; use BlockMixing for $n groups"))

    M = zeros(Float64, 3, 3)
    M[1, 1] = r[1] * activity[1]
    M[1, 2] = (1 - r[1]) * activity[1] * split[1]
    M[1, 3] = (1 - r[1]) * activity[1] * (1 - split[1])

    M[2, 1] = (1 - r[2]) * activity[2] * (1 - split[2])
    M[2, 2] = r[2] * activity[2]
    M[2, 3] = (1 - r[2]) * activity[2] * split[2]

    M[3, 1] = (1 - r[3]) * activity[3] * (1 - split[3])
    M[3, 2] = (1 - r[3]) * activity[3] * split[3]
    M[3, 3] = r[3] * activity[3]
    M
end

function _lift_block(diagonal::AbstractMatrix{<:Real}, offdiagonal::AbstractMatrix{<:Real},
                     target_base::Int, source_base::Int)
    if target_base == source_base
        diagonal
    elseif target_base > source_base
        offdiagonal
    else
        transpose(offdiagonal)
    end
end
