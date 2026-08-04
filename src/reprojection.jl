"""
Demographic reprojection morphism.

Adapts a contact matrix built for one population to a different population at the
*same* partition — the population-changing counterpart to coarsening/refinement,
which change partition at fixed population. The formula is the reciprocity-preserving
projection of Arregui, Aleta, Sanz & Moreno, "Projecting social contact matrices to
different demographic structures", PLoS Comput Biol 14(12):e1006638 (2018), and
matches `socialmixr::contact_matrix(..., symmetric = TRUE)`'s
`normalise_weighted_matrix` when its `survey_pop` argument differs from the survey's
own population.
"""

"""
    reproject(cm::ContactMatrix, target_population::AbstractVector{<:Real})

Reproject `cm` onto `target_population`, the same groups under a different
demographic structure.

Written in `MeanContacts`, where the formula is

    M'[i,j] = (M[i,j]·N'_j + M[j,i]·N'_i) / (2·N'_j)

i.e. `symmetrise` evaluated at `N'` instead of `cm`'s own population — the source
population plays no role in the formula, only in getting `cm` into `MeanContacts`
first when it isn't already there. Carried to the other representations by the
representation isomorphisms, so that

    reinterpret_units(reproject(cm, N'), s) == reproject(reinterpret_units(cm, s), N')

Reciprocity preservation is unconditional in every representation: `reproject(cm,
N')` is always reciprocal under `N'`, regardless of whether `cm` was reciprocal
under its own population.

`reproject` is **not functorial** in the population
(`reproject(reproject(cm, N'), N'') != reproject(cm, N'')` in general — reprojecting
through an intermediate population differs from reprojecting directly to the final
one) and it **does not commute** with `coarsen` (reprojecting then coarsening
differs from coarsening then reprojecting). 

If a group has zero target population, a reciprocal finite rate only exists when
the corresponding total contacts are also zero; otherwise an `ArgumentError` is
thrown, matching `symmetrise`.

# Returns
A new `ContactMatrix` in the same representation as `cm`, carrying
`target_population` and reciprocal under it.
"""
reproject(cm::ContactMatrix, target_population::AbstractVector{<:Real}) =
    _via(MeanContacts(), c -> _reproject_mean(c, target_population), cm)

# Reprojection in its home representation, `MeanContacts`. Call `reproject` rather
# than this: applying it directly to another representation is the error `_via`
# exists to prevent.
function _reproject_mean(cm::ContactMatrix, target_population::AbstractVector{<:Real})
    n = n_groups(cm)
    length(target_population) == n || throw(DimensionMismatch(
        "target population length $(length(target_population)) does not match partition with $n groups"))
    N′ = Float64.(collect(target_population))
    M_proj = _symmetrise_at(matrix(cm), N′)
    ContactMatrix(M_proj, cm.partition, N′, cm.semantics)
end
