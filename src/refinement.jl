"""
Refinement: splitting a coarse contact matrix over a finer partition.

Unlike coarsening (which is canonical), refinement needs auxiliary assumptions —
here a population prior giving each fine group's share of its coarse group.

Refinement is a section of coarsening, not an inverse: coarsening a refined matrix
returns the original, while refining a coarsened one does not, since coarsening
discards the within-group structure the prior then has to supply.
"""

"""
    refine(cm::ContactMatrix, fine::AbstractPartition, fine_population::AbstractVector{<:Real})

Refine a coarse contact matrix to a finer partition, splitting each coarse group's
contacts between its fine subgroups in proportion to their populations:

    M_fine[i, j] = M_coarse[I, J] * N_i / N_I

for fine groups `i, j` lying in coarse groups `I, J`, where `N_I` is the total
population of coarse group `I`.

Only the contactee (row) index carries a population share. In total contacts the
split is proportional on both axes; the participant-side share cancels against the
division by participant population that turns total contacts into mean contacts.

Coarsening the result along the same map recovers `cm`'s matrix (proved in
`proofs/ContACTProofs/Refinement.lean`). Total contacts are preserved, and the
result is reciprocal whenever `cm` is — exactly for a prior that aggregates exactly,
and otherwise up to the tolerance on that check.

A coarse group with zero population supplies no shares to split by, so it is
accepted only when it carries no contacts on either axis. A fine group with zero
population comes out zero on both axes — its row by the formula, its column by
construction.

# Arguments
- `cm`: coarse contact matrix
- `fine`: target fine partition (must be a refinement of `cm.partition`)
- `fine_population`: population per fine group, which must aggregate to `cm`'s
  population within `rtol=1e-8, atol=1e-8`, and must agree exactly on which coarse
  groups are empty. This is the prior that determines the split

Written in `MeanContacts`, where the split is the single factor `N_i / N_I`, and
carried to the other representations by the representation isomorphisms, so that
`reinterpret_units(refine(cm, fine, p), s) ≈ refine(reinterpret_units(cm, s), fine, p)`
— exact in real arithmetic, agreeing to about 15 significant digits in floating point.
"""
refine(cm::ContactMatrix, fine::AbstractPartition,
       fine_population::AbstractVector{<:Real}) =
    _via(MeanContacts(), c -> _refine_mean(c, fine, fine_population), cm)

# Refinement in its home representation, `MeanContacts`.
function _refine_mean(cm::ContactMatrix, fine::AbstractPartition,
                      fine_population::AbstractVector{<:Real})
    coarse = cm.partition
    n_fine = n_groups(fine)
    n_coarse = n_groups(coarse)

    length(fine_population) == n_fine || throw(DimensionMismatch(
        "fine_population length $(length(fine_population)) ≠ fine groups $n_fine"))

    # Build the mapping from fine to coarse
    f = PartitionMap(fine, coarse)
    fmap = _assignments(f)

    M_coarse = matrix(cm)
    pop_fine = Float64.(fine_population)

    # Population within each coarse group
    pop_coarse = zeros(Float64, n_coarse)
    for i in 1:n_fine
        pop_coarse[fmap[i]] += pop_fine[i]
    end
    pop_cm = population(cm)
    isapprox(pop_coarse, pop_cm; rtol=1e-8, atol=1e-8) || throw(ArgumentError(
        "fine_population must aggregate to the ContactMatrix population under the refinement map"))

    # That check scales its tolerance by the largest population, so it would let a
    # prior put whole people into an empty coarse group — bypassing the guard below
    # on the one group it exists to protect. Which groups are empty must match exactly.
    for I in 1:n_coarse
        (pop_coarse[I] == 0) == (pop_cm[I] == 0) || throw(ArgumentError(
            "coarse group $I is empty in one population but not the other " *
            "(matrix $(pop_cm[I]), aggregated fine_population $(pop_coarse[I])); " *
            "fine_population must sum to exactly zero for an empty coarse group"))
    end

    for I in 1:n_coarse
        pop_coarse[I] == 0 || continue
        all(iszero, @view M_coarse[I, :]) || throw(ArgumentError(
            "cannot refine coarse group $I: zero population with nonzero contacts in its row, " *
            "so there are no population shares to split them by"))
        all(iszero, @view M_coarse[:, I]) || throw(ArgumentError(
            "cannot refine coarse group $I: zero population with nonzero contacts in its column, " *
            "which its fine subgroups cannot carry"))
    end

    # Each fine group's share of its coarse group's contacts. The guard above leaves
    # only zero rows at a zero coarse population, so the 0.0 branch multiplies zeros.
    share = [pop_coarse[fmap[i]] > 0 ? pop_fine[i] / pop_coarse[fmap[i]] : 0.0
             for i in 1:n_fine]

    M_fine = zeros(Float64, n_fine, n_fine)
    for j in 1:n_fine
        pop_fine[j] == 0 && continue   # no participants, so the column stays zero
        J = fmap[j]
        for i in 1:n_fine
            M_fine[i, j] = M_coarse[fmap[i], J] * share[i]
        end
    end

    ContactMatrix(M_fine, fine, pop_fine, cm.semantics)
end
