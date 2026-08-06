"""
Unit semantics as a groupoid of representations.

The `UnitSemantics` tags do not name different objects. They name *representations*
of one object, related by invertible row and column rescalings. The quantity they
all scale is the **total contacts** `T[i,j] = M[i,j]·N_j`, the number of contacts
group-`j` participants report with group-`i` contactees. A representation is a pair
of exponents `(a, b)`:

    R[i,j] = T[i,j] · N_i^(-a) · N_j^(-b)

| Semantics       | `(a,b)` | entry `R[i,j]` | reading                          |
|-----------------|---------|----------------|----------------------------------|
| `ContactCounts` | `(0,0)` | `M[i,j]·N_j`   | total contacts between the pair  |
| `MeanContacts`  | `(0,1)` | `M[i,j]`       | per participant                  |
| `PerCapitaRate` | `(1,1)` | `M[i,j]/N_i`   | per participant, per contactee   |

The exponents are constrained, not free: a per-capita denominator has to range over
something a contact event actually involves, and one contact involves exactly one
participant (in `j`) and one contactee (in `i`). So the admissible representations
are the corners of the unit square `a, b ∈ {0,1}`. `M[i,j]/N_j` — dividing by the
participant population a second time — is outside it, and is not a rate of anything.

`PerCapitaRate` matches `socialmixr::per_capita`, which divides by the **contactee**
population. socialmixr indexes rows by participant and columns by contact, the
transpose of ContACT, so its "divide each column by that group's population" is
ContACT's *row* division. Copying that rule without transposing is what produced the
earlier `M[i,j]/N_j` definition; the same trap as the Arregui transpose caveat.

## Reciprocity depends only on `d = a - b`

Reciprocity is symmetry of total contacts, so in representation `(a,b)` it reads

    R[i,j] · N_i^d == R[j,i] · N_j^d,        d = a - b

which is `R[i,j]·N_j == R[j,i]·N_i` for `MeanContacts` (`d = -1`) and plain matrix
symmetry for both `ContactCounts` and `PerCapitaRate` (`d = 0`). The row and column
directions cancel in `d`, so the reciprocity law cannot see the component of a
rescaling along `a = b` — see the naturality theorems in
`proofs/ContACTProofs/UnitSemantics.lean`.

That `PerCapitaRate` is symmetric under reciprocity is the property the object is
named for: Lomas et al. (2025)'s `β_ij = a_i·a_j/D` and Wallinga et al. (2006)'s
`β_ij = q·m_ij/w_i` are both this representation.

**Domain restriction.** A row or column with zero population can only be rescaled
when it is already zero. Conversion throws otherwise, in either direction, so that
the conversions really are mutually inverse on their common domain.
"""

"""
    population_exponents(s::UnitSemantics)

Exponents `(a, b)` placing `s` in the representation lattice: an entry in
representation `s` is `T[i,j]·N_i^(-a)·N_j^(-b)`, where `T[i,j] = M[i,j]·N_j` is
total contacts.
"""
population_exponents(::ContactCounts) = (0, 0)
population_exponents(::MeanContacts) = (0, 1)
population_exponents(::PerCapitaRate) = (1, 1)

"""
    reinterpret_units(cm::ContactMatrix, target::UnitSemantics)

Convert `cm` to the `target` representation, dispatching on the semantics `cm`
actually carries.

With `(a,b)` the exponents of `cm.semantics` and `(a′,b′)` those of `target`, row
`i` is multiplied by `N_i^(a-a′)` and column `j` by `N_j^(b-b′)`. Converting to the
representation a matrix is already in is the identity.

Throws an `ArgumentError` if a group with zero population carries nonzero contacts
in a direction that is being rescaled, since no rescaling of that row or column is
invertible.

See also [`to_counts`](@ref), [`to_mean_contacts`](@ref), [`to_per_capita`](@ref).
"""
function reinterpret_units(cm::ContactMatrix, target::UnitSemantics)
    a, b = population_exponents(cm.semantics)
    a′, b′ = population_exponents(target)
    row_exp = a - a′
    col_exp = b - b′

    M = matrix(cm)
    pop = population(cm)
    n = n_groups(cm)

    (row_exp == 0 && col_exp == 0) &&
        return ContactMatrix(M, cm.partition, pop, target)

    if row_exp != 0
        for i in 1:n
            pop[i] == 0 && !all(iszero, @view M[i, :]) && throw(ArgumentError(
                "cannot reinterpret row $i: zero population with nonzero contacts"))
        end
    end
    if col_exp != 0
        for j in 1:n
            pop[j] == 0 && !all(iszero, @view M[:, j]) && throw(ArgumentError(
                "cannot reinterpret column $j: zero population with nonzero contacts"))
        end
    end

    R = zeros(Float64, n, n)
    for j in 1:n
        (pop[j] == 0 && col_exp != 0) && continue   # column is zero; leave it zero
        for i in 1:n
            (pop[i] == 0 && row_exp != 0) && continue
            R[i, j] = M[i, j] * Float64(pop[i])^row_exp * Float64(pop[j])^col_exp
        end
    end
    ContactMatrix(R, cm.partition, pop, target)
end

# Transport `f` into the representation `cm` is in.
#
# A structural morphism is written once, in the representation where it is simplest —
# its *home* — and defined everywhere else by conjugation with the representation
# isomorphisms. This is what makes every morphism natural with respect to
# `reinterpret_units`:
#
#     reinterpret_units(f(cm), s) ≈ f(reinterpret_units(cm, s))
#
# The law is an equality of contact matrices, but it is written with `≈`, and the
# docstrings of the morphisms that rely on it follow suit, for two reasons. The two
# routes multiply and divide by the populations in a different order, so their
# entries agree to about 15 significant digits rather than bit for bit — close
# enough that `==` on the numbers is false more often than true, but well inside
# `≈`'s default relative tolerance of `sqrt(eps())`, which is roughly 8 significant
# digits. A discrepancy anywhere near that tolerance is a bug, not rounding. And
# `ContactMatrix` defines no `==` of its own, so Julia compares the wrapped arrays
# by identity: two matrices with identical contents are already `!=` before any
# arithmetic enters into it.
#
# `f` receives a `ContactMatrix` tagged `home` and must return one tagged `home`;
# carrying `cm.semantics` forward, as the implementations do, achieves that. The
# result is converted back using the result's own population, so morphisms that
# change partition and population (`↓`, `↑`, `⊗`) transport correctly.
#
# When `cm` is already in `home` this is exactly `f(cm)`.
function _via(home::UnitSemantics, f, cm::ContactMatrix)
    cm.semantics isa typeof(home) && return f(cm)
    reinterpret_units(f(reinterpret_units(cm, home)), cm.semantics)
end

"""
    to_mean_contacts(cm::ContactMatrix)

Convert to `MeanContacts`: the mean number of group-`i` contacts reported by a
group-`j` participant. This is the representation nearly every operation in the
package is stated in.
"""
to_mean_contacts(cm::ContactMatrix) = reinterpret_units(cm, MeanContacts())

"""
    to_per_capita(cm::ContactMatrix)

Convert to `PerCapitaRate`, dividing row `i` by the population of the contactee
(row) group so that entries are per participant *and* per contactee.

This matches `socialmixr::per_capita` once the row/column conventions are lined up,
and is the object the reciprocity literature calls a per-capita contact rate or a
transmission rate `β_ij`: Wallinga et al. (2006)'s `β_ij = q·m_ij/w_i` and Lomas et
al. (2025)'s `β_ij = a_i·a_j/D` are both this representation. For a reciprocal
matrix the result is **symmetric**, which is the property the object is named for.

Dispatches on the semantics of `cm`, so converting from `ContactCounts` divides by
both populations rather than assuming `MeanContacts` input.
"""
to_per_capita(cm::ContactMatrix) = reinterpret_units(cm, PerCapitaRate())

"""
    to_counts(cm::ContactMatrix)

Convert to `ContactCounts`, the total contacts reported between each ordered pair
of groups, multiplying column `j` by the population of the participant (column)
group.
"""
to_counts(cm::ContactMatrix) = reinterpret_units(cm, ContactCounts())
