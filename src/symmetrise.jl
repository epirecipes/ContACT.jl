"""
Symmetrisation morphism.

Makes a contact matrix reciprocal, so that the total contacts group `i` reports
with group `j` equal the total group `j` reports with group `i`.

Key properties (proven in Lean): symmetrisation establishes reciprocity, is
idempotent, and is **natural** with respect to the unit-semantics isomorphisms.
Naturality is not re-derived per representation — the morphism is defined once in
its home representation and conjugated by the isomorphisms.
"""

"""
    symmetrise(cm::ContactMatrix)

Symmetrise a contact matrix, making total contacts reciprocal.

Written in `MeanContacts`, where the formula is

    M_sym[i, j] = (M[i,j]·N_j + M[j,i]·N_i) / (2·N_j)

and carried to the other representations by the representation isomorphisms, so that

    reinterpret_units(symmetrise(cm), s) == symmetrise(reinterpret_units(cm, s))

for every representation `s`. Reciprocity itself is representation-dependent —
`M[i,j]·N_j == M[j,i]·N_i` under `MeanContacts`, plain matrix symmetry under
`ContactCounts` and `PerCapitaRate` — and transport is what reconciles the two
facts: one formula, the right answer everywhere.

If a group has zero population, reciprocal finite rates only exist when the
corresponding total contacts are also zero; otherwise an `ArgumentError` is thrown.

# Returns
A new `ContactMatrix` in the same representation, with a reciprocal contact pattern.
"""
symmetrise(cm::ContactMatrix) = _via(MeanContacts(), _symmetrise_mean, cm)

# Symmetrisation in its home representation, `MeanContacts`. Call `symmetrise`
# rather than this: applying it directly to another representation is the error
# `_via` exists to prevent.
_symmetrise_mean(cm::ContactMatrix) =
    ContactMatrix(_symmetrise_at(matrix(cm), population(cm)), cm.partition, population(cm), cm.semantics)

# Force an n×n `MeanContacts`-style matrix reciprocal against `pop`.
function _symmetrise_at(M::AbstractMatrix{<:Real}, pop::AbstractVector{<:Real})
    n = size(M, 1)
    M_sym = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            # Total contacts from i→j and j→i, averaged
            total_ij = M[i, j] * pop[j] + M[j, i] * pop[i]
            if pop[j] == 0
                total_ij == 0 || throw(ArgumentError(
                    "cannot symmetrise pair ($i, $j): group $j has zero population but nonzero total contacts"))
                M_sym[i, j] = 0.0
            else
                M_sym[i, j] = total_ij / (2.0 * pop[j])
            end
        end
    end
    M_sym
end
