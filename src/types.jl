"""
Core types for ContACT.jl.

The fundamental objects in the categorical framework:
- `AgePartition`: describes how a continuous age range is discretized
- `ContactSurvey`: individual-level contact data (participants + contacts)
- `ContactMatrix`: the categorical object combining matrix + partition + population + semantics
"""

# ---------------------------------------------------------------------------
# Unit semantics: what the matrix entries represent
# ---------------------------------------------------------------------------

"""
    UnitSemantics

Abstract type for the semantic interpretation of contact matrix entries.
"""
abstract type UnitSemantics end

"""Mean number of contacts reported by participants in row group with contacts in column group."""
struct MeanContacts <: UnitSemantics end

"""Raw count of contacts between age groups (not normalised by participants)."""
struct ContactCounts <: UnitSemantics end

"""Per-capita contact rate: mean contacts divided by population of contact group."""
struct PerCapitaRate <: UnitSemantics end

# ---------------------------------------------------------------------------
# Partitions
# ---------------------------------------------------------------------------

"""
    AbstractPartition{D}

Finite partition of survey records into groups along semantic dimension `D`.

The dimension is a symbol (for example `:age`, `:sex`, or `:region`) or a tuple
of symbols for product partitions. Concrete partitions know which participant
and contact columns to read from a `ContactSurvey`.
"""
abstract type AbstractPartition{D} end

"""Semantic dimension of a partition."""
dimension(::AbstractPartition{D}) where {D} = D

_default_participant_col(d::Symbol) = Symbol(:part_, d)
_default_contact_col(d::Symbol) = Symbol(:cnt_, d)

"""
    IntervalPartition{D,T}

Partition a numeric variable into intervals defined by lower limits.

`AgePartition` is an alias for `IntervalPartition{:age,Float64}`.
"""
struct IntervalPartition{D,T<:Real} <: AbstractPartition{D}
    participant_col::Symbol
    contact_col::Symbol
    limits::Vector{T}
    labels::Vector{String}
end

function IntervalPartition{D,T}(limits::AbstractVector{<:Real};
                                participant_col::Symbol=_default_participant_col(D),
                                contact_col::Symbol=_default_contact_col(D),
                                labels::Vector{String}=String[]) where {D,T<:Real}
    raw_limits = collect(T, limits)
    sorted = sort(raw_limits)
    length(sorted) > 0 || throw(ArgumentError("partition limits must be non-empty"))
    all(isfinite, sorted) || throw(ArgumentError("partition limits must be finite"))
    all(>=(zero(T)), sorted) || throw(ArgumentError("partition limits must be non-negative"))
    allunique(sorted) || throw(ArgumentError("partition limits must be unique"))

    if isempty(labels)
        labels = _make_interval_labels(sorted)
    else
        length(labels) == length(sorted) ||
            throw(ArgumentError("length(labels) must equal length(limits)"))
        raw_limits == sorted ||
            throw(ArgumentError("partition limits must already be sorted when custom labels are supplied"))
    end
    IntervalPartition{D,T}(participant_col, contact_col, sorted, labels)
end

IntervalPartition{D}(limits::AbstractVector{<:Real}; kwargs...) where {D} =
    IntervalPartition{D,Float64}(limits; kwargs...)

IntervalPartition(dimension::Symbol, limits::AbstractVector{<:Real}; kwargs...) =
    IntervalPartition{dimension,Float64}(limits; kwargs...)

IntervalPartition(dimension::Symbol; breaks::AbstractVector{<:Real}, kwargs...) =
    IntervalPartition(dimension, breaks; kwargs...)

"""
    AgePartition

Compatibility alias for `IntervalPartition{:age,Float64}` using participant
column `:part_age` and contact column `:cnt_age` by default.
"""
const AgePartition = IntervalPartition{:age,Float64}

"""
    CategoricalPartition{D,T}

Partition a categorical variable into a fixed set of levels.
"""
struct CategoricalPartition{D,T} <: AbstractPartition{D}
    participant_col::Symbol
    contact_col::Symbol
    levels::Vector{T}
    labels::Vector{String}
    lookup::Dict{T,Int}
end

function CategoricalPartition{D,T}(levels::AbstractVector;
                                   participant_col::Symbol=_default_participant_col(D),
                                   contact_col::Symbol=_default_contact_col(D),
                                   labels::Vector{String}=String[]) where {D,T}
    levs = collect(T, levels)
    length(levs) > 0 || throw(ArgumentError("categorical levels must be non-empty"))
    any(ismissing, levs) && throw(ArgumentError("categorical levels must not contain missing"))
    allunique(levs) || throw(ArgumentError("categorical levels must be unique"))
    if isempty(labels)
        labels = string.(levs)
    else
        length(labels) == length(levs) ||
            throw(ArgumentError("length(labels) must equal length(levels)"))
    end
    lookup = Dict(level => i for (i, level) in enumerate(levs))
    CategoricalPartition{D,T}(participant_col, contact_col, levs, labels, lookup)
end

CategoricalPartition{D}(levels::AbstractVector; kwargs...) where {D} =
    CategoricalPartition{D,eltype(levels)}(levels; kwargs...)

CategoricalPartition(dimension::Symbol, levels::AbstractVector; kwargs...) =
    CategoricalPartition{dimension,eltype(levels)}(levels; kwargs...)

CategoricalPartition(dimension::Symbol; levels::AbstractVector, kwargs...) =
    CategoricalPartition(dimension, levels; kwargs...)

"""
    ProductPartition(factors...)

Cartesian product of finite partitions, e.g. `age × sex`.
"""
struct ProductPartition{Ds,P<:Tuple} <: AbstractPartition{Ds}
    factors::P
    labels::Vector{String}
end

function ProductPartition(factors::Tuple)
    length(factors) > 0 || throw(ArgumentError("product partition needs at least one factor"))
    all(f -> f isa AbstractPartition, factors) ||
        throw(ArgumentError("all product factors must be partitions"))
    dims = Tuple(dimension(f) for f in factors)
    allunique(dims) || throw(ArgumentError(
        "product partition dimensions must be unique; got $dims"))
    labels = _product_labels(Tuple(group_labels(f) for f in factors))
    ProductPartition{dims,typeof(factors)}(factors, labels)
end

ProductPartition(factors::AbstractPartition...) = ProductPartition(factors)

×(a::AbstractPartition, b::AbstractPartition) = ProductPartition(a, b)
×(a::ProductPartition, b::AbstractPartition) = ProductPartition(a.factors..., b)
×(a::AbstractPartition, b::ProductPartition) = ProductPartition(a, b.factors...)
×(a::ProductPartition, b::ProductPartition) = ProductPartition(a.factors..., b.factors...)

"""Number of groups in the partition."""
n_groups(p::AbstractPartition) = length(group_labels(p))

"""Human-readable group labels."""
group_labels(p::AbstractPartition) = p.labels

"""Compatibility alias for existing age-oriented API."""
age_labels(p::AbstractPartition) = group_labels(p)

"""Lower interval limits. Defined for interval partitions."""
age_limits(p::IntervalPartition) = p.limits
age_limits(p::AbstractPartition) =
    throw(ArgumentError("age_limits is only defined for interval partitions; use group_labels for general partitions"))

function _make_interval_labels(limits::Vector{<:Real})
    n = length(limits)
    labels = String[]
    for i in 1:n
        lo = Int(limits[i])
        if i < n
            hi = Int(limits[i + 1])
            push!(labels, "[$lo,$hi)")
        else
            push!(labels, "$lo+")
        end
    end
    return labels
end

function _product_labels(label_sets::Tuple)
    isempty(label_sets) && return [""]
    first_labels = first(label_sets)
    rest_labels = _product_labels(Base.tail(label_sets))
    return [isempty(rest) ? first : string(first, ":", rest)
            for first in first_labels for rest in rest_labels]
end

function assign_group(p::IntervalPartition, value)
    value === missing && return nothing
    x = Float64(value)
    isfinite(x) || return nothing
    grp = searchsortedlast(p.limits, x)
    grp == 0 && return nothing
    return grp
end

function assign_group(p::CategoricalPartition, value)
    value === missing && return nothing
    return get(p.lookup, value, nothing)
end

function _cartesian_index(indices::Tuple, sizes::Tuple)
    idx = 1
    stride = 1
    for k in length(indices):-1:1
        idx += (indices[k] - 1) * stride
        stride *= sizes[k]
    end
    return idx
end

function _cartesian_indices(index::Integer, sizes::Tuple)
    1 <= index <= prod(sizes) || throw(BoundsError(1:prod(sizes), index))
    rem = index - 1
    indices = Vector{Int}(undef, length(sizes))
    for k in length(sizes):-1:1
        indices[k] = rem % sizes[k] + 1
        rem ÷= sizes[k]
    end
    Tuple(indices)
end

function assign_participant_group(p::Union{IntervalPartition,CategoricalPartition}, row)
    hasproperty(row, p.participant_col) ||
        throw(ArgumentError("participants must have :$(p.participant_col) column for $(dimension(p)) partition"))
    assign_group(p, row[p.participant_col])
end

function assign_contact_group(p::Union{IntervalPartition,CategoricalPartition}, row)
    hasproperty(row, p.contact_col) ||
        throw(ArgumentError("contacts must have :$(p.contact_col) column for $(dimension(p)) partition"))
    assign_group(p, row[p.contact_col])
end

function assign_participant_group(p::ProductPartition, row)
    groups = Tuple(assign_participant_group(f, row) for f in p.factors)
    any(g -> g === nothing, groups) && return nothing
    _cartesian_index(groups, Tuple(n_groups(f) for f in p.factors))
end

function assign_contact_group(p::ProductPartition, row)
    groups = Tuple(assign_contact_group(f, row) for f in p.factors)
    any(g -> g === nothing, groups) && return nothing
    _cartesian_index(groups, Tuple(n_groups(f) for f in p.factors))
end

same_partition(a::IntervalPartition{D}, b::IntervalPartition{D}) where {D} =
    a.limits == b.limits && a.labels == b.labels

same_partition(a::CategoricalPartition{D}, b::CategoricalPartition{D}) where {D} =
    a.levels == b.levels && a.labels == b.labels

same_partition(a::ProductPartition{D}, b::ProductPartition{D}) where {D} =
    length(a.factors) == length(b.factors) &&
    all(same_partition(x, y) for (x, y) in zip(a.factors, b.factors))

same_partition(::AbstractPartition, ::AbstractPartition) = false

Base.:(==)(a::AbstractPartition, b::AbstractPartition) = same_partition(a, b)

# ---------------------------------------------------------------------------
# Contact survey
# ---------------------------------------------------------------------------

"""
    ContactSurvey

Individual-level social contact data: participants who reported contacts.

# Fields
- `participants`: DataFrame with at least column `:part_id`
- `contacts`: DataFrame with at least column `:part_id`
- `metadata`: Dict of additional information (country, year, reference, etc.)
"""
struct ContactSurvey
    participants::DataFrame
    contacts::DataFrame
    metadata::Dict{Symbol, Any}

    function ContactSurvey(participants::DataFrame, contacts::DataFrame;
                           metadata::Dict{Symbol, Any}=Dict{Symbol, Any}())
        :part_id in propertynames(participants) ||
            throw(ArgumentError("participants must have :part_id column"))
        :part_id in propertynames(contacts) ||
            throw(ArgumentError("contacts must have :part_id column"))
        _validate_contact_survey_ids(participants, contacts)
        new(participants, contacts, metadata)
    end
end

function _validate_contact_survey_ids(participants::DataFrame, contacts::DataFrame)
    participant_ids = participants[!, :part_id]
    any(ismissing, participant_ids) &&
        throw(ArgumentError("participants.part_id must not contain missing values"))
    allunique(participant_ids) ||
        throw(ArgumentError("participants.part_id values must be unique"))

    contact_ids = contacts[!, :part_id]
    any(ismissing, contact_ids) &&
        throw(ArgumentError("contacts.part_id must not contain missing values"))
    participant_id_set = Set(participant_ids)
    unknown = count(id -> !(id in participant_id_set), contact_ids)
    unknown == 0 || throw(ArgumentError(
        "contacts contain $unknown part_id value(s) not present in participants"))
    nothing
end

# ---------------------------------------------------------------------------
# Contact matrix — the main categorical object
# ---------------------------------------------------------------------------

"""
    ContactMatrix{T<:Real, S<:UnitSemantics}

A contact matrix bundled with its finite partition, population vector, and
unit semantics. This is the object type in the category **Contact**.

The matrix entry `M[i,j]` represents interactions from group `j`
(column = "contactor") to group `i` (row = "contacted"), with
interpretation given by `S`.

# Fields
- `matrix`: n×n contact matrix
- `partition`: finite group partition
- `population`: population size per group
- `semantics`: what entries represent (MeanContacts, ContactCounts, PerCapitaRate)
"""
struct ContactMatrix{T<:Real, S<:UnitSemantics, P<:AbstractPartition}
    matrix::Matrix{T}
    partition::P
    population::Vector{T}
    semantics::S

    function ContactMatrix(M::AbstractMatrix{<:Real}, partition::P,
                           population::AbstractVector{<:Real},
                           semantics::S=MeanContacts()) where {P<:AbstractPartition, S<:UnitSemantics}
        n = n_groups(partition)
        size(M) == (n, n) || throw(DimensionMismatch(
            "matrix size $(size(M)) does not match partition with $n groups"))
        length(population) == n || throw(DimensionMismatch(
            "population length $(length(population)) does not match partition with $n groups"))
        all(x -> isfinite(x) && x >= 0, M) ||
            throw(ArgumentError("contact matrix entries must be finite and non-negative"))
        all(x -> isfinite(x) && x >= 0, population) ||
            throw(ArgumentError("population entries must be finite and non-negative"))

        T = promote_type(eltype(M), eltype(population))
        pop = Vector{T}(population)
        new{T, S, P}(Matrix{T}(M), partition, pop, semantics)
    end
end

"""Extract the raw matrix."""
matrix(cm::ContactMatrix) = cm.matrix

"""Extract the population vector."""
population(cm::ContactMatrix) = cm.population

"""Number of partition groups."""
n_groups(cm::ContactMatrix) = n_groups(cm.partition)

"""Partition labels and compatibility accessors."""
age_limits(cm::ContactMatrix) = age_limits(cm.partition)
age_labels(cm::ContactMatrix) = age_labels(cm.partition)
group_labels(cm::ContactMatrix) = group_labels(cm.partition)

function Base.show(io::IO, cm::ContactMatrix{T, S}) where {T, S}
    n = n_groups(cm)
    print(io, "ContactMatrix{$T, $(nameof(S))} ($n×$n groups)")
end
