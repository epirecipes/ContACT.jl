"""
Coarsening: pushforward of contact matrices along finite partition maps.

Categorically, given a surjective map f: Fine → Coarse between finite group
sets, coarsening is the left Kan extension along f. It preserves total contacts.

The key property (proven in Lean): coarsen(g ∘ f) = coarsen(g) ∘ coarsen(f)
"""

"""
    PartitionMap

A morphism between partitions of the same semantic dimension. It maps each group
in the domain partition to a group in the codomain partition and wraps a Catlab
`FinFunction`.

`AgeMap` is a compatibility alias for `PartitionMap{:age}`.

# Example
```julia
fine = AgePartition([0, 5, 10, 15, 20, 65])
coarse = AgePartition([0, 15, 65])
f = PartitionMap(fine, coarse, [1, 1, 1, 2, 2, 3])
```
"""
struct PartitionMap{D,E,P<:AbstractPartition{D},Q<:AbstractPartition{E}}
    domain::P
    codomain::Q
    mapping::FinFunction
end

function _partition_map(domain::P, codomain::Q,
                        assignments::AbstractVector{Int}) where {D,E,P<:AbstractPartition{D},Q<:AbstractPartition{E}}
    n_domain = n_groups(domain)
    n_codomain = n_groups(codomain)
    length(assignments) == n_domain || throw(ArgumentError(
        "assignments length $(length(assignments)) ≠ domain groups $n_domain"))
    all(1 .<= assignments .<= n_codomain) || throw(ArgumentError(
        "all assignments must be in 1:$n_codomain"))
    PartitionMap{D,E,P,Q}(domain, codomain, FinFunction(assignments, n_codomain))
end

function PartitionMap{D,E}(domain::P, codomain::Q,
                           assignments::AbstractVector{Int}) where {D,E,P<:AbstractPartition{D},Q<:AbstractPartition{E}}
    _partition_map(domain, codomain, assignments)
end

PartitionMap(domain::AbstractPartition{D}, codomain::AbstractPartition{E},
             assignments::AbstractVector{Int}) where {D,E} =
    _partition_map(domain, codomain, assignments)

"""
    AgeMap

Compatibility alias for `PartitionMap{:age}`.
"""
const AgeMap = PartitionMap{:age,:age}

AgeMap(domain::AgePartition, codomain::AgePartition) = PartitionMap(domain, codomain)
AgeMap(domain::AgePartition, codomain::AgePartition, assignments::AbstractVector{Int}) =
    _partition_map(domain, codomain, assignments)

"""
    PartitionMap(fine::IntervalPartition, coarse::IntervalPartition)

Automatically construct a partition map from compatible interval partitions
where coarse limits are a subset of fine limits.
"""
function PartitionMap{D}(fine::IntervalPartition{D}, coarse::IntervalPartition{D}) where {D}
    fine_lims = age_limits(fine)
    coarse_lims = age_limits(coarse)
    all(l ∈ fine_lims for l in coarse_lims) || throw(ArgumentError(
        "coarse interval limits must be a subset of fine interval limits"))

    assignments = Int[]
    coarse_idx = 1
    for (i, l) in enumerate(fine_lims)
        if coarse_idx < length(coarse_lims) && i > 1 && l >= coarse_lims[coarse_idx + 1]
            coarse_idx += 1
        end
        push!(assignments, coarse_idx)
    end
    PartitionMap{D,D}(fine, coarse, assignments)
end

PartitionMap(fine::IntervalPartition{D}, coarse::IntervalPartition{D}) where {D} =
    PartitionMap{D}(fine, coarse)

"""
    PartitionMap(domain::CategoricalPartition, codomain::CategoricalPartition, map)

Construct a categorical coarsening map from a dictionary mapping domain levels to
codomain levels.
"""
function PartitionMap{D}(domain::CategoricalPartition{D}, codomain::CategoricalPartition{D},
                         level_map::AbstractDict) where {D}
    assignments = Int[]
    for level in domain.levels
        haskey(level_map, level) ||
            throw(ArgumentError("missing mapping for level $(repr(level))"))
        target = level_map[level]
        idx = assign_group(codomain, target)
        idx === nothing &&
            throw(ArgumentError("mapped level $(repr(target)) is not in codomain"))
        push!(assignments, idx)
    end
    PartitionMap{D,D}(domain, codomain, assignments)
end

PartitionMap(domain::CategoricalPartition{D}, codomain::CategoricalPartition{D},
             level_map::AbstractDict) where {D} =
    PartitionMap{D}(domain, codomain, level_map)

function PartitionMap{D}(domain::CategoricalPartition{D}, codomain::CategoricalPartition{D}) where {D}
    assignments = Int[]
    for level in domain.levels
        idx = assign_group(codomain, level)
        idx === nothing &&
            throw(ArgumentError("cannot infer categorical map: level $(repr(level)) is not in codomain"))
        push!(assignments, idx)
    end
    PartitionMap{D,D}(domain, codomain, assignments)
end

PartitionMap(domain::CategoricalPartition{D}, codomain::CategoricalPartition{D}) where {D} =
    PartitionMap{D}(domain, codomain)

"""
    PartitionMap(product::ProductPartition, factor::AbstractPartition)

Construct the projection from a product partition to one of its factors.
"""
function PartitionMap(domain::ProductPartition, codomain::AbstractPartition)
    matches = findall(f -> same_partition(f, codomain), collect(domain.factors))
    length(matches) == 1 || throw(ArgumentError(
        "cannot infer product projection: codomain must match exactly one product factor"))

    factor_idx = only(matches)
    sizes = Tuple(n_groups(f) for f in domain.factors)
    assignments = [_cartesian_indices(i, sizes)[factor_idx] for i in 1:n_groups(domain)]
    PartitionMap(domain, codomain, assignments)
end

"""
    coarsen(cm::ContactMatrix, f::PartitionMap)

Coarsen a contact matrix along a finite partition map (left Kan extension).

The coarsened matrix preserves total contacts: for each coarse pair (I, J),
the entry is the population-weighted average of fine entries mapping to (I, J).

Specifically:
    M_coarse[I, J] = Σ_{j∈f⁻¹(J)} (N_j / N_J) * Σ_{i∈f⁻¹(I)} M_fine[i, j]

where N_j is the population of fine group j and N_J = Σ_{j∈f⁻¹(J)} N_j.
"""
function coarsen(cm::ContactMatrix, f::PartitionMap)
    same_partition(cm.partition, f.domain) || throw(ArgumentError(
        "ContactMatrix partition does not match PartitionMap domain"))

    M = matrix(cm)
    pop = population(cm)
    n_fine = n_groups(cm)
    n_coarse = n_groups(f.codomain)
    fmap = collect(f.mapping)

    pop_coarse = zeros(Float64, n_coarse)
    for j in 1:n_fine
        pop_coarse[fmap[j]] += pop[j]
    end

    M_coarse = zeros(Float64, n_coarse, n_coarse)
    for j in 1:n_fine
        J = fmap[j]
        weight_j = pop_coarse[J] > 0 ? pop[j] / pop_coarse[J] : 0.0
        for i in 1:n_fine
            I = fmap[i]
            M_coarse[I, J] += weight_j * M[i, j]
        end
    end

    ContactMatrix(M_coarse, f.codomain, pop_coarse, cm.semantics)
end

"""
    coarsen(cm::ContactMatrix, coarse::AbstractPartition)

Coarsen a contact matrix to a coarser partition by auto-constructing a
`PartitionMap` when possible.
"""
function coarsen(cm::ContactMatrix, coarse::AbstractPartition)
    f = PartitionMap(cm.partition, coarse)
    coarsen(cm, f)
end
