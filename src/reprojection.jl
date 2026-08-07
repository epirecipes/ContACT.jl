"""
Demographic reprojection morphism.

Adapts a contact matrix built for one population to a different population at the
*same* partition — the population-changing counterpart to coarsening/refinement,
which change partition at fixed population. The formula averages each pair's two
reported contact directions under the new population, from Arregui, Aleta, Sanz &
Moreno, "Projecting social contact matrices to different demographic structures",
PLoS Comput Biol 14(12):e1006638 (2018), and matches
`socialmixr::contact_matrix(..., symmetric = TRUE)`'s `normalise_weighted_matrix`
when its `survey_pop` argument differs from the survey's own population.

That paper compares four projection methods. This one restores reciprocity but
changes the mixing pattern while doing so, by an amount that grows with how far the
target demography is from the survey's; see `reproject`'s docstring.
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

    reinterpret_units(reproject(cm, N'), s) ≈ reproject(reinterpret_units(cm, s), N')

Reciprocity is established unconditionally in every representation: `reproject(cm,
N')` is always reciprocal under `N'`, whether or not `cm` was reciprocal under its
own population. On input that already is, reprojection preserves that reciprocity
across the population change; on input that is not, it repairs it.

`reproject` is **not functorial** in the population: reprojecting through an
intermediate population differs from reprojecting straight to the final one, so
`reproject(reproject(cm, N'), N'')` and `reproject(cm, N'')` disagree entrywise in
general. It also **does not commute** with `coarsen` — reprojecting then coarsening
differs from coarsening then reprojecting.

If a group has zero target population, a reciprocal finite rate only exists when
the corresponding total contacts are also zero; otherwise an `ArgumentError` is
thrown, matching `symmetrise`.

# Accuracy

Averaging the two directions restores reciprocity, but it does not correct contact
densities for the change in demographic structure, and Arregui et al. quantify the
cost. Projecting a 2005 Polish survey onto 2018 population data, they report the
ratio of projected to original contact density exceeding 1.5 for some age pairs and
falling to about 0.5 for others — ">50% bias in the contact densities projected
between certain age groups" — from ignoring how the age-strata fractions shifted
over those 13 years. In their SEIR forecasts this method also predicts incidence
trends that their density-corrected alternative substantially reduces or reverses.

Leaving the matrix uncorrected is worse. Hamilton, Knight & Mishra (2024) show that
imbalanced matrices always underestimate R₀, measuring 3.1–5.7% bias in their
most-imbalanced countries and up to ~39% bias in downstream vaccination impact,
using this same correction. So the question these figures bear on is how far a
matrix can be projected, not whether to correct it.

ContACT does not implement the density-corrected alternatives, as right now there is no
rule for how a density-corrected matrix should combine under `↓`, `↑` and `⊗`.

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
