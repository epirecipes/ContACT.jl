"""
Source-stratified rectangular contact matrices for partial survey observations.

These matrices support workflows such as CoMix/SEP, where participant-side
covariates are richer than contact-side covariates. Rows are contacted/target
groups and columns are contactor/source groups, matching `ContactMatrix`
orientation.
"""

"""
    SourceStratifiedContactMatrix

Rectangular contact matrix whose columns are source/participant groups and rows
are target/contact groups.

This is useful when contact-side covariates are only partially observed. For
example, a survey may identify participant groups by `age × SEP × education`
while contacts are reported only by age. Entry `M[i, j]` is the mean number of
contacts in target group `i` reported by participants in source group `j`.
"""
struct SourceStratifiedContactMatrix{T<:Real,S<:UnitSemantics,
                                     R<:AbstractPartition,C<:AbstractPartition}
    matrix::Matrix{T}
    target_partition::R
    source_partition::C
    source_population::Vector{T}
    semantics::S

    function SourceStratifiedContactMatrix(M::AbstractMatrix{<:Real},
                                           target_partition::R,
                                           source_partition::C,
                                           source_population::AbstractVector{<:Real},
                                           semantics::S=MeanContacts()) where {
                                               R<:AbstractPartition,
                                               C<:AbstractPartition,
                                               S<:UnitSemantics,
                                           }
        n_target = n_groups(target_partition)
        n_source = n_groups(source_partition)
        size(M) == (n_target, n_source) || throw(DimensionMismatch(
            "matrix size $(size(M)) does not match target/source groups ($n_target, $n_source)"))
        length(source_population) == n_source || throw(DimensionMismatch(
            "source_population length $(length(source_population)) does not match source groups $n_source"))
        all(x -> isfinite(x) && x >= 0, M) ||
            throw(ArgumentError("source-stratified matrix entries must be finite and non-negative"))
        all(x -> isfinite(x) && x >= 0, source_population) ||
            throw(ArgumentError("source population entries must be finite and non-negative"))

        T = promote_type(eltype(M), eltype(source_population))
        new{T,S,R,C}(
            Matrix{T}(M),
            target_partition,
            source_partition,
            Vector{T}(source_population),
            semantics,
        )
    end
end

"""Extract the raw rectangular matrix."""
matrix(cm::SourceStratifiedContactMatrix) = cm.matrix

"""Extract the source/participant population vector."""
population(cm::SourceStratifiedContactMatrix) = cm.source_population

"""Extract the target/contact partition."""
target_partition(cm::SourceStratifiedContactMatrix) = cm.target_partition

"""Extract the source/participant partition."""
source_partition(cm::SourceStratifiedContactMatrix) = cm.source_partition

"""Number of target/contact groups."""
n_target_groups(cm::SourceStratifiedContactMatrix) = n_groups(cm.target_partition)

"""Number of source/participant groups."""
n_source_groups(cm::SourceStratifiedContactMatrix) = n_groups(cm.source_partition)

"""Human-readable target/contact group labels."""
target_group_labels(cm::SourceStratifiedContactMatrix) = group_labels(cm.target_partition)

"""Human-readable source/participant group labels."""
source_group_labels(cm::SourceStratifiedContactMatrix) = group_labels(cm.source_partition)

Base.size(cm::SourceStratifiedContactMatrix) = size(matrix(cm))
Base.size(cm::SourceStratifiedContactMatrix, d::Integer) = size(matrix(cm), d)

function Base.show(io::IO, cm::SourceStratifiedContactMatrix{T,S}) where {T,S}
    print(io, "SourceStratifiedContactMatrix{$T, $(nameof(S))} ",
        "($(n_target_groups(cm))×$(n_source_groups(cm)) target×source groups)")
end

"""
    compute_source_stratified_matrix(survey, target_partition, source_partition;
                                     population=nothing, weights=nothing)

Compute a rectangular source-stratified contact matrix from partial survey data.

The `source_partition` is applied to participant rows and the `target_partition`
is applied to contact rows. This allows participant/source columns to include
richer covariates than are available for contacts, e.g. columns indexed by
`age × SEP × education` and rows indexed by age only.

Mean contacts are always normalised by the number of *surveyed* participants in
each source group. The stored `source_population`, however, is the supplied
`population` vector when given (otherwise the surveyed counts). When a census
`population` is passed, `source_total_contacts` therefore reports post-stratified
totals (mean contacts × census population), which differ from the raw observed
totals; a source group with surveyed participants but zero supplied population
contributes zero total contacts.
"""
function compute_source_stratified_matrix(
    survey::ContactSurvey,
    target_partition::AbstractPartition,
    source_partition::AbstractPartition;
    population::Union{Nothing,AbstractVector{<:Real}}=nothing,
    weights::Union{Nothing,Symbol}=nothing,
)
    n_target = n_groups(target_partition)
    n_source = n_groups(source_partition)

    source_groups = [assign_participant_group(source_partition, row)
                     for row in eachrow(survey.participants)]
    target_groups = [assign_contact_group(target_partition, row)
                     for row in eachrow(survey.contacts)]
    w = _survey_participant_weights(survey, weights)

    contact_counts = zeros(Float64, n_target, n_source)
    participant_weights = zeros(Float64, n_source)
    part_id_to_idx = Dict(id => i for (i, id) in enumerate(survey.participants.part_id))
    dropped_missing_contact_group = 0
    dropped_unknown_participant = 0
    dropped_missing_participant_group = 0

    for (k, target_group) in enumerate(target_groups)
        if target_group === nothing
            dropped_missing_contact_group += 1
            continue
        end
        participant_idx = get(part_id_to_idx, survey.contacts.part_id[k], nothing)
        if participant_idx === nothing
            dropped_unknown_participant += 1
            continue
        end
        source_group = source_groups[participant_idx]
        if source_group === nothing
            dropped_missing_participant_group += 1
            continue
        end
        contact_counts[target_group, source_group] += w[participant_idx]
    end
    _warn_dropped_contacts("compute_source_stratified_matrix";
        missing_contact_group=dropped_missing_contact_group,
        unknown_participant=dropped_unknown_participant,
        missing_participant_group=dropped_missing_participant_group)

    for (i, source_group) in enumerate(source_groups)
        source_group === nothing && continue
        participant_weights[source_group] += w[i]
    end

    M = zeros(Float64, n_target, n_source)
    for source_group in 1:n_source
        if participant_weights[source_group] > 0
            M[:, source_group] .= contact_counts[:, source_group] ./ participant_weights[source_group]
        end
    end

    source_pop = if population === nothing
        participant_weights
    else
        Float64.(population)
    end

    SourceStratifiedContactMatrix(
        M,
        target_partition,
        source_partition,
        source_pop,
        MeanContacts(),
    )
end

function _survey_participant_weights(survey::ContactSurvey, weights::Union{Nothing,Symbol})
    if weights === nothing
        return ones(Float64, nrow(survey.participants))
    end

    weights in propertynames(survey.participants) ||
        throw(ArgumentError("weight column :$weights not found in participants"))
    raw_weights = survey.participants[!, weights]
    any(ismissing, raw_weights) &&
        throw(ArgumentError("weight column :$weights contains missing values"))
    w = Float64.(raw_weights)
    all(x -> isfinite(x) && x >= 0, w) ||
        throw(ArgumentError("weight column :$weights must be finite and non-negative"))
    w
end

"""
    source_total_contacts(cm::SourceStratifiedContactMatrix)

Return total contacts by multiplying each source column by its source
population.
"""
function source_total_contacts(cm::SourceStratifiedContactMatrix)
    M = matrix(cm)
    pop = population(cm)
    counts = zeros(Float64, size(M))
    for source in 1:n_source_groups(cm)
        counts[:, source] .= M[:, source] .* pop[source]
    end
    counts
end

"""
    coarsen_sources(cm::SourceStratifiedContactMatrix, f::PartitionMap)

Coarsen the source/participant columns of a rectangular matrix while preserving
total contacts. The target/contact rows are unchanged.
"""
function coarsen_sources(cm::SourceStratifiedContactMatrix, f::PartitionMap)
    same_partition(cm.source_partition, f.domain) || throw(ArgumentError(
        "SourceStratifiedContactMatrix source partition does not match PartitionMap domain"))

    n_target = n_target_groups(cm)
    n_source_coarse = n_groups(f.codomain)
    fmap = _assignments(f)
    source_pop = population(cm)

    coarse_pop = zeros(Float64, n_source_coarse)
    for source in 1:n_source_groups(cm)
        coarse_pop[fmap[source]] += source_pop[source]
    end

    counts = source_total_contacts(cm)
    coarse_counts = zeros(Float64, n_target, n_source_coarse)
    for source in 1:n_source_groups(cm)
        coarse_counts[:, fmap[source]] .+= counts[:, source]
    end

    coarse_matrix = _source_mean_contacts_from_counts(coarse_counts, coarse_pop)
    SourceStratifiedContactMatrix(
        coarse_matrix,
        cm.target_partition,
        f.codomain,
        coarse_pop,
        cm.semantics,
    )
end

# Rescale `fine_pop` so that each fibre of `fmap` sums to its base group's population.
# Each fibre is compared with its own base group, so a small group is held to the same
# relative tolerance as a large one, and the two must agree exactly on which groups are
# empty. A fibre that already sums exactly is returned unchanged.
function _fibre_matched_population(fine_pop::AbstractVector{<:Real},
                                   fmap::AbstractVector{Int},
                                   base_pop::AbstractVector{<:Real},
                                   name::AbstractString)
    aggregated = zeros(Float64, length(base_pop))
    for (k, I) in enumerate(fmap)
        aggregated[I] += fine_pop[k]
    end
    for I in eachindex(base_pop)
        (aggregated[I] == 0) == (base_pop[I] == 0) || throw(ArgumentError(
            "base group $I is empty in one population but not the other " *
            "(base population $(base_pop[I]), $name summed over the group $(aggregated[I]))"))
        isapprox(aggregated[I], base_pop[I]; rtol=1e-8) || throw(ArgumentError(
            "$name must sum to the base population in every base group; group $I sums to " *
            "$(aggregated[I]) against a base population of $(base_pop[I])"))
    end

    matched = Vector{Float64}(fine_pop)
    for (k, I) in enumerate(fmap)
        aggregated[I] == base_pop[I] && continue
        matched[k] *= base_pop[I] / aggregated[I]
    end
    matched
end

# Total contacts `M[i,j]·N_j` of a base matrix that a lift splits into blocks, returned
# exactly symmetric. Each pair `(i, j)` is checked at its own scale, so a small group's
# contacts are held to the same relative tolerance as a large group's. A zero-population
# group must have a zero participant column, since its column carries no total contacts
# and cannot be recovered by coarsening the lift.
function _lift_base_counts(cm::ContactMatrix, operation::AbstractString)
    M = matrix(cm)
    pop = population(cm)
    n = n_groups(cm)
    for j in 1:n
        pop[j] == 0 || continue
        all(iszero, @view M[:, j]) || throw(ArgumentError(
            "$operation cannot use base group $j: zero population with nonzero contacts " *
            "in its participant column"))
    end

    C = M * Diagonal(pop)
    for j in 1:n, i in 1:(j - 1)
        a, b = C[i, j], C[j, i]
        abs(a - b) <= 1e-8 * max(abs(a), abs(b)) || throw(ArgumentError(
            "$operation requires reciprocal total contacts (M[i,j]·N_j == M[j,i]·N_i), " *
            "but groups $i and $j give $a and $b. If rows may be participants rather " *
            "than contactees, transpose the matrix; otherwise apply ↔ first"))
    end
    (C .+ transpose(C)) ./ 2
end

"""
    align_source_stratified_matrix(partial, base[, f])

Scale a source-stratified intermediate matrix so that coarsening its source
groups matches a reciprocal base `ContactMatrix` in total contacts and in mean
contacts.

The map `f` must send `partial.source_partition` to `base.partition`. If omitted,
ContACT attempts to infer it, for example as the projection
`age × SEP × education → age`.

# Requirements
- `base` must be reciprocal in total contacts, `M[i,j]·N_j ≈ M[j,i]·N_i` to a
  relative tolerance of `1e-8` for every pair. Apply `↔(base)` first when
  starting from a raw survey estimate. A zero-population base group must have a
  zero participant column.
- The source population of `partial` must sum, over the source groups in each
  base group, to that group's base population to a relative tolerance of `1e-8`,
  and must sum to exactly zero for an empty base group.

# Returns
A `SourceStratifiedContactMatrix` whose total contacts, summed over the source
groups in each base group, equal the symmetrised base total contacts
`(C + Cᵀ)/2`, where `C[i,j] = M[i,j]·N_j`. A base total of at most `1e-8` for
which `partial` observes no contacts is left at zero; a larger one throws. The
source population is `partial`'s, rescaled within each base group to sum to the
base population.
"""
function align_source_stratified_matrix(
    partial::SourceStratifiedContactMatrix,
    base::ContactMatrix,
    f::PartitionMap=PartitionMap(partial.source_partition, base.partition),
)
    partial.semantics isa MeanContacts || throw(ArgumentError(
        "source-stratified alignment requires MeanContacts semantics"))
    base.semantics isa MeanContacts || throw(ArgumentError(
        "source-stratified alignment requires a MeanContacts base matrix"))
    same_partition(partial.target_partition, base.partition) || throw(ArgumentError(
        "partial target partition must match the base ContactMatrix partition"))
    same_partition(partial.source_partition, f.domain) || throw(ArgumentError(
        "source map domain must match the partial source partition"))
    same_partition(base.partition, f.codomain) || throw(ArgumentError(
        "source map codomain must match the base ContactMatrix partition"))

    fmap = _assignments(f)
    n_target = n_target_groups(partial)
    n_base = n_groups(base.partition)
    source_pop = population(partial)
    base_pop = population(base)

    aligned_pop = _fibre_matched_population(source_pop, fmap, base_pop, "source population")
    base_counts = _lift_base_counts(base, "source-stratified alignment")

    counts = source_total_contacts(partial)
    current = zeros(Float64, n_target, n_base)
    for source in 1:n_source_groups(partial)
        current[:, fmap[source]] .+= counts[:, source]
    end

    aligned_counts = copy(counts)
    for target in 1:n_target
        for base_source in 1:n_base
            desired = base_counts[target, base_source]
            observed = current[target, base_source]
            if observed > 0
                scale = desired / observed
                for source in 1:n_source_groups(partial)
                    fmap[source] == base_source || continue
                    aligned_counts[target, source] *= scale
                end
            elseif desired > 1e-8
                throw(ArgumentError(
                    "cannot align target group $target and source base group $base_source: observed total is zero but base total is $desired"))
            end
        end
    end

    aligned_matrix = _source_mean_contacts_from_counts(aligned_counts, aligned_pop)
    SourceStratifiedContactMatrix(
        aligned_matrix,
        partial.target_partition,
        partial.source_partition,
        aligned_pop,
        MeanContacts(),
    )
end

function _source_mean_contacts_from_counts(counts::AbstractMatrix{<:Real},
                                           pop::AbstractVector{<:Real})
    M = zeros(Float64, size(counts))
    for source in eachindex(pop)
        if pop[source] == 0
            all(iszero, counts[:, source]) || throw(ArgumentError(
                "cannot convert source counts for zero-population source group $source"))
        else
            M[:, source] .= counts[:, source] ./ pop[source]
        end
    end
    M
end
