"""
Coarsening: pushforward of contact matrices along age-group aggregation maps.

Categorically: given a surjective map f: Fine → Coarse between age partitions,
coarsening is the left Kan extension along f. It preserves total contacts.

The key property (proven in Lean): coarsen(g ∘ f) = coarsen(g) ∘ coarsen(f)
"""

"""
    AgeMap

A morphism between age partitions: maps each fine age group to a coarse group.
Wraps a Catlab `FinFunction`.

# Example
```julia
fine = AgePartition([0, 5, 10, 15, 20, 65])
coarse = AgePartition([0, 15, 65])
# Groups 1,2,3 → 1; Group 4 → 2; Group 5 → 2; Group 6 → 3
f = AgeMap(fine, coarse, [1, 1, 1, 2, 2, 3])
```
"""
struct AgeMap
    domain::AgePartition
    codomain::AgePartition
    mapping::FinFunction

    function AgeMap(domain::AgePartition, codomain::AgePartition,
                    assignments::AbstractVector{Int})
        n_fine = n_groups(domain)
        n_coarse = n_groups(codomain)
        length(assignments) == n_fine || throw(ArgumentError(
            "assignments length $(length(assignments)) ≠ domain groups $n_fine"))
        all(1 .<= assignments .<= n_coarse) || throw(ArgumentError(
            "all assignments must be in 1:$n_coarse"))
        f = FinFunction(assignments, n_coarse)
        new(domain, codomain, f)
    end
end

"""
    AgeMap(fine::AgePartition, coarse::AgePartition)

Automatically construct an AgeMap from compatible partitions where coarse
limits are a subset of fine limits.
"""
function AgeMap(fine::AgePartition, coarse::AgePartition)
    fine_lims = age_limits(fine)
    coarse_lims = age_limits(coarse)
    all(l ∈ fine_lims for l in coarse_lims) || throw(ArgumentError(
        "coarse limits must be a subset of fine limits"))

    assignments = Int[]
    coarse_idx = 1
    for (i, l) in enumerate(fine_lims)
        if coarse_idx < length(coarse_lims) && i > 1 && l >= coarse_lims[coarse_idx + 1]
            coarse_idx += 1
        end
        push!(assignments, coarse_idx)
    end
    AgeMap(fine, coarse, assignments)
end

"""
    coarsen(cm::ContactMatrix, f::AgeMap)

Coarsen a contact matrix along an age-group map (left Kan extension).

The coarsened matrix preserves total contacts: for each coarse pair (I, J),
the entry is the population-weighted average of fine entries mapping to (I, J).

Specifically:
    M_coarse[I, J] = Σ_{j∈f⁻¹(J)} (N_j / N_J) * Σ_{i∈f⁻¹(I)} M_fine[i, j]

where N_j is the population of fine group j and N_J = Σ_{j∈f⁻¹(J)} N_j.
"""
function coarsen(cm::ContactMatrix, f::AgeMap)
    cm.partition.limits == f.domain.limits || throw(ArgumentError(
        "ContactMatrix partition does not match AgeMap domain"))

    M = matrix(cm)
    pop = population(cm)
    n_fine = n_groups(cm)
    n_coarse = n_groups(f.codomain)
    fmap = collect(f.mapping)

    # Compute coarse population
    pop_coarse = zeros(Float64, n_coarse)
    for j in 1:n_fine
        pop_coarse[fmap[j]] += pop[j]
    end

    # Compute coarsened matrix
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
    coarsen(cm::ContactMatrix, coarse::AgePartition)

Coarsen a contact matrix to a coarser partition (auto-constructs the AgeMap).
"""
function coarsen(cm::ContactMatrix, coarse::AgePartition)
    f = AgeMap(cm.partition, coarse)
    coarsen(cm, f)
end
