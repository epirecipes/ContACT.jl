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
# Age partition
# ---------------------------------------------------------------------------

"""
    AgePartition

Discretisation of age into groups defined by lower limits.

# Fields
- `limits`: sorted vector of lower age limits (e.g., [0, 5, 15, 65])
- `labels`: human-readable labels (auto-generated if not provided)
"""
struct AgePartition
    limits::Vector{Float64}
    labels::Vector{String}

    function AgePartition(limits::AbstractVector{<:Real}; labels::Vector{String}=String[])
        sorted = sort(collect(Float64, limits))
        length(sorted) > 0 || throw(ArgumentError("age limits must be non-empty"))
        allunique(sorted) || throw(ArgumentError("age limits must be unique"))

        if isempty(labels)
            labels = _make_age_labels(sorted)
        else
            length(labels) == length(sorted) ||
                throw(ArgumentError("length(labels) must equal length(limits)"))
        end
        new(sorted, labels)
    end
end

"""Number of age groups in the partition."""
n_groups(p::AgePartition) = length(p.limits)

"""Lower age limits."""
age_limits(p::AgePartition) = p.limits

"""Human-readable labels."""
age_labels(p::AgePartition) = p.labels

function _make_age_labels(limits::Vector{Float64})
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

# ---------------------------------------------------------------------------
# Contact survey
# ---------------------------------------------------------------------------

"""
    ContactSurvey

Individual-level social contact data: participants who reported contacts.

# Fields
- `participants`: DataFrame with at least columns `:part_id` and `:part_age`
- `contacts`: DataFrame with at least columns `:part_id` and `:cnt_age`
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
        :part_age in propertynames(participants) ||
            throw(ArgumentError("participants must have :part_age column"))
        :part_id in propertynames(contacts) ||
            throw(ArgumentError("contacts must have :part_id column"))
        :cnt_age in propertynames(contacts) ||
            throw(ArgumentError("contacts must have :cnt_age column"))
        new(participants, contacts, metadata)
    end
end

# ---------------------------------------------------------------------------
# Contact matrix — the main categorical object
# ---------------------------------------------------------------------------

"""
    ContactMatrix{T<:Real, S<:UnitSemantics}

A contact matrix bundled with its age partition, population vector, and
unit semantics. This is the object type in the category **Contact**.

The matrix entry `M[i,j]` represents interactions from age group `j`
(column = "contactor") to age group `i` (row = "contacted"), with
interpretation given by `S`.

# Fields
- `matrix`: n×n contact matrix
- `partition`: age group discretisation
- `population`: population size per age group
- `semantics`: what entries represent (MeanContacts, ContactCounts, PerCapitaRate)
"""
struct ContactMatrix{T<:Real, S<:UnitSemantics}
    matrix::Matrix{T}
    partition::AgePartition
    population::Vector{T}
    semantics::S

    function ContactMatrix(M::AbstractMatrix{T}, partition::AgePartition,
                           population::AbstractVector{<:Real},
                           semantics::S=MeanContacts()) where {T<:Real, S<:UnitSemantics}
        n = n_groups(partition)
        size(M) == (n, n) || throw(DimensionMismatch(
            "matrix size $(size(M)) does not match partition with $n groups"))
        length(population) == n || throw(DimensionMismatch(
            "population length $(length(population)) does not match partition with $n groups"))
        pop = Vector{T}(population)
        new{T, S}(Matrix{T}(M), partition, pop, semantics)
    end
end

"""Extract the raw matrix."""
matrix(cm::ContactMatrix) = cm.matrix

"""Extract the population vector."""
population(cm::ContactMatrix) = cm.population

"""Number of age groups."""
n_groups(cm::ContactMatrix) = n_groups(cm.partition)

"""Age partition."""
age_limits(cm::ContactMatrix) = age_limits(cm.partition)
age_labels(cm::ContactMatrix) = age_labels(cm.partition)

function Base.show(io::IO, cm::ContactMatrix{T, S}) where {T, S}
    n = n_groups(cm)
    print(io, "ContactMatrix{$T, $(nameof(S))} ($n×$n age groups)")
end
