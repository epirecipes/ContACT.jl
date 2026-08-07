"""
Population transport morphism.

Changes a contact matrix's population at a *fixed* partition while holding total
contacts fixed. `reproject` also changes the population at a fixed partition, but
the two preserve different things, and only `transport_population` is
path-independent under successive population changes. Both are perfectly composable
as functions; what `reproject` fails is the population-change composition law,
`R_{N'→N''} ∘ R_{N→N'} == R_{N→N''}`.

`transport_population(cm, N')` keeps `cm`'s own total contacts `M[i,j]·N_j` and
changes only the population that turns them back into a rate. `reproject(cm, N')`
instead reads `cm`'s entries as though they had been measured under `N'`, then
averages the two directions to make them reciprocal; the source population never
enters its formula. So transport carries the source matrix's contacts forward
unchanged, while reprojection recomputes them under the new population and repairs
the result.

That is also why transport satisfies identity and composition laws that
reprojection cannot. Repairing reciprocity is not the identity on a matrix that
was not already reciprocal, so `reproject` cannot be the identity when the
population does not change.
"""

"""
    transport_population(cm::ContactMatrix, target_population::AbstractVector{<:Real})

Transport `cm` onto `target_population`, the same groups under a different
population, holding total contacts fixed.

Written in `MeanContacts`, where the formula is

    M'[i,j] = M[i,j]·N_j/N'_j

and carried to the other representations by the representation isomorphisms. On the
common domain of those conversions — in particular whenever populations are strictly
positive —

    reinterpret_units(transport_population(cm, N'), s) ≈
        transport_population(reinterpret_units(cm, s), N')

At a zero population the conversions have different domains, so this can fail; see
the note on the zero-source guard below.

The quantity held fixed is total contacts, so in `ContactCounts` this operation is
the identity on the matrix and changes only the attached population. It is
nevertheless implemented in `MeanContacts`, where nearly every other operation in
the package lives, so the common path requires no unit conversion.

Because coarsening only sums total contacts over the fibres of a partition map,
`transport_population` commutes with coarsening whenever the coarse target
population is obtained by aggregating the fine target population along the same map.
This holds for any valid fibre map; the fine groups need not be contiguous or
order-preserving. Proved as `weightedCoarsen_populationTransport_nonneg` in
`proofs/ContACTProofs/PopulationTransport.lean` — the variant whose hypotheses match
this function's domain, since `ContactMatrix` populations are nonnegative by
construction. (`weightedCoarsen_populationTransport` is the same statement over
unrestricted reals, where it needs every coarse population nonzero.) Checked against
this implementation by the `"transport: coarsening equivariance"` testset.

Four laws hold, stated here on matrix entries:

    identity      M'[i,j] == M[i,j]                    when N' == N
    composition   transporting N → N' → N'' agrees with N → N'' directly
    inverse       transporting N → N' → N recovers M
    totals        M'[i,j]·N'_j == M[i,j]·N_j

All four are exact in real arithmetic. In floating point the last three agree to
about 15 significant digits, so tests assert them with `≈`. The identity law is
bitwise exact: a group whose population does not change has its column copied rather
than rescaled, since `M·N/N` can round away from `M`. Transporting onto `cm`'s own
population returns an equal matrix, provided that population is strictly positive —
the positivity check runs first, so it throws otherwise. The result is a fresh
`ContactMatrix`, not `cm` itself, so mutating it cannot reach back into `cm`.

If `cm` is already reciprocal at its own population, the result is reciprocal at
`target_population`. Unlike `reproject`, this operation does not repair
non-reciprocal input — it carries any imbalance through unchanged, in
total-contacts terms.

`target_population` must be strictly positive: a zero target makes the division
singular, and an empty target group has no well-defined rate even when its total
contacts vanish, since transport does not touch total contacts at all. For an
empty target group, `reproject` is available when its zero-total condition is
satisfied.

A zero *source* population is admitted only when its participant column is zero, so
that transport never silently discards contacts from an empty participant group.
Such a transport is one-way — it cannot be undone, since reversing it would need a
zero target — but it still composes with later transports onto positive populations,
and the composition law holds across it. The identity and inverse laws, and hence the
groupoid reading, are what require strictly positive populations throughout.

That guard is on the column alone, so it is *not* the same domain as the
representation isomorphisms, which guard whichever axis they rescale. Converting
between `MeanContacts` and `PerCapitaRate` rescales rows, so a matrix whose empty
group has a zero column but a nonzero row is accepted here under `MeanContacts` or
`ContactCounts` and rejected under `PerCapitaRate`. Which matrices are accepted at
a zero population therefore does depend on the representation.

# Returns
A new `ContactMatrix` in the same representation as `cm`, carrying
`target_population`. Always a new object, never `cm` itself.
"""
function transport_population(cm::ContactMatrix, target_population::AbstractVector{<:Real})
    n = n_groups(cm)
    length(target_population) == n || throw(DimensionMismatch(
        "target population length $(length(target_population)) does not match partition with $n groups"))
    N′ = Float64.(collect(target_population))
    all(>(0), N′) || throw(ArgumentError(
        "transport_population requires a strictly positive target population " *
        "(got a non-positive entry); use reproject for target populations with " *
        "empty groups"))

    # An unchanged population scales nothing, so short-circuit rather than send the
    # matrix through a representation round trip, which would round in each
    # direction. Rebuild rather than return `cm`: every operation here hands back
    # arrays independent of its input.
    N′ == population(cm) &&
        return ContactMatrix(matrix(cm), cm.partition, population(cm), cm.semantics)

    _via(MeanContacts(), c -> _transport_population_mean(c, N′), cm)
end

# Population transport in its home representation, `MeanContacts`. Call
# `transport_population` rather than this: applying it directly to another
# representation is the error `_via` exists to prevent, and the target population
# arrives already validated.
_transport_population_mean(cm::ContactMatrix, target_population::Vector{Float64}) =
    ContactMatrix(_transport_at(matrix(cm), population(cm), target_population),
                  cm.partition, target_population, cm.semantics)

# Scale each column `j` of an n×n `MeanContacts`-style matrix by
# `source_pop[j]/target_pop[j]`, which holds that column's total contacts fixed.
function _transport_at(M::AbstractMatrix{<:Real}, source_pop::AbstractVector{<:Real},
                        target_pop::AbstractVector{<:Real})
    n = size(M, 1)

    # An empty participant group has no contacts to carry, so a nonzero column is
    # rejected rather than silently zeroed. The guard is on the column only, which
    # is narrower than what `reinterpret_units` checks: that guards whichever axis
    # it rescales, so a zero-population group with a zero column but a nonzero row
    # passes here and still fails a conversion to `PerCapitaRate`.
    for j in 1:n
        source_pop[j] == 0 && !all(iszero, @view M[:, j]) && throw(ArgumentError(
            "cannot transport column $j: zero population with nonzero contacts"))
    end

    M_t = zeros(Float64, n, n)
    for j in 1:n
        if source_pop[j] == target_pop[j]
            # Copy rather than multiply and divide by the same number, which is
            # not the identity in floating point.
            for i in 1:n
                M_t[i, j] = M[i, j]
            end
        else
            for i in 1:n
                M_t[i, j] = M[i, j] * source_pop[j] / target_pop[j]
            end
        end
    end
    M_t
end
