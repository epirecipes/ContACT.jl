"""
Refinement: parameterized disaggregation of contact matrices.

Unlike coarsening (which is canonical), refinement requires auxiliary
assumptions — typically a population prior for proportional splitting.
This is NOT claimed to be an inverse of coarsening.
"""

"""
    refine(cm::ContactMatrix, fine::AbstractPartition, fine_population::AbstractVector{<:Real})

Refine a coarse contact matrix to a finer partition using proportional
disaggregation based on population distribution within each coarse group.

Each coarse entry M[I, J] is distributed to fine entries proportionally:
    M_fine[i, j] = M_coarse[I, J]  (contacts per participant unchanged)

The fine population determines how participants are distributed among
sub-groups within each coarse bin.

# Arguments
- `cm`: coarse contact matrix
- `fine`: target fine partition (must be a refinement of `cm.partition`)
- `fine_population`: population per fine age group (the distributional prior)
"""
function refine(cm::ContactMatrix, fine::AbstractPartition,
                fine_population::AbstractVector{<:Real})
    coarse = cm.partition
    n_fine = n_groups(fine)
    n_coarse = n_groups(coarse)

    length(fine_population) == n_fine || throw(DimensionMismatch(
        "fine_population length $(length(fine_population)) ≠ fine groups $n_fine"))

    # Build the mapping from fine to coarse
    f = PartitionMap(fine, coarse)
    fmap = collect(f.mapping)

    M_coarse = matrix(cm)
    pop_fine = Float64.(fine_population)

    # Population within each coarse group
    pop_coarse = zeros(Float64, n_coarse)
    for i in 1:n_fine
        pop_coarse[fmap[i]] += pop_fine[i]
    end
    isapprox(pop_coarse, population(cm); rtol=1e-8, atol=1e-8) || throw(ArgumentError(
        "fine_population must aggregate to the ContactMatrix population under the refinement map"))

    # Proportional disaggregation
    M_fine = zeros(Float64, n_fine, n_fine)
    for j in 1:n_fine
        J = fmap[j]
        for i in 1:n_fine
            I = fmap[i]
            # Contact rate is assumed uniform within coarse bins
            M_fine[i, j] = M_coarse[I, J]
        end
    end

    ContactMatrix(M_fine, fine, pop_fine, cm.semantics)
end
