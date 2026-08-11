# ContACT.jl

**Category theory applied to contact matrices** — a Julia package for building, composing, and manipulating structured contact matrices using applied category theory.

## Overview

ContACT.jl provides a categorical framework for:
- Computing contact matrices from individual-level survey data over arbitrary finite partitions (▷)
- Composing setting-specific matrices (home, work, school) via **additive composition** (⊕)
- Changing partition resolution via **coarsening** (population-weighted pushforward) and **refinement** (its section, given a population prior)
- Refining reciprocal matrices by respondent social activity, following Britton-Ball-style high/low or quantile activity strata (⤊)
- Lifting age matrices to generalized age × SES/contact matrices with explicit mixing assumptions (⊠)
- Reconstructing generalized matrices from partial participant-side covariates with source-stratified intermediate matrices
- Adding spatial structure via **stratification** (Kronecker product) (⊗)
- Enforcing reciprocity via **symmetrisation** (↔)

Structural operations are verified with formal Lean 4 proofs of their claimed algebraic laws (see the proof table below), while numerical utilities, such as spectral radius, R₀ helpers, and epidemic bounds, undergo standard testing.

## Quick Start

```julia
using ContACT
using DataFrames

# Define an age partition (age is one interval-valued partition)
partition = AgePartition([0, 5, 18, 65])

# Or define a categorical partition using symbols from the survey data
sex = CategoricalPartition(:sex;
    participant_col=:part_sex,
    contact_col=:cnt_sex,
    levels=["F", "M"],
)

# Compute contact matrix from survey data
cm = compute_matrix(survey, partition; population=pop)

# Compose setting-specific matrices
total = cm_home ⊕ cm_work ⊕ cm_school ⊕ cm_other

# Coarsen to fewer age groups
coarse = AgePartition([0, 18, 65])
cm_coarse = cm ↓ coarse

# Stratify across 3 regions
coupling = [0.8 0.1 0.1; 0.1 0.8 0.1; 0.1 0.1 0.8]
cm_spatial = cm ⊗ coupling

# Symmetrise (enforce reciprocity)
cm_sym = ↔(cm)

# Lift to age × activity strata using respondent contact rates
spec = ActivityRefinement(survey; n=2, mixing=:proportionate)
cm_activity = cm_sym ⤊ spec

# Lift to age × SES strata using an explicit generalized-matrix assumption
ses = CategoricalPartition(:ses; levels=["low", "middle", "high"])
gcm_spec = GeneralizedLift(ses; distribution=[0.35, 0.45, 0.20])
cm_age_ses = cm_sym ⊠ gcm_spec

# For partial data, build an intermediate matrix with richer participant columns
partial = compute_source_stratified_matrix(survey, partition, partition × ses)
```

## Operators

| Symbol | LaTeX | Operation | Structure |
|--------|-------|-----------|-------------------|
| `⊕` | `\oplus` | Additive composition | Commutative monoid at fixed partition, population, and semantics |
| `⊗` | `\otimes` | Stratification | Kronecker product with a coupling matrix |
| `↓` | `\downarrow` | Coarsening | Population-weighted pushforward |
| `↑` | `\uparrow` | Population-based refinement | Section of coarsening |
| `⤊` | `\Uuparrow` | Activity refinement | Hidden-stratum lift |
| `⊠` | `\boxtimes` | Generalized product lift | Assumption-driven product lift |
| `▷` | `\triangleright` | Survey-to-matrix construction | Normalisation of survey contact counts by participant mass |
| `∘` | `\circ` | PartitionMap composition | Partition-map composition |
| `↔` | `\leftrightarrow` | Symmetrisation | Reciprocity projection |
| `ρ` | `\rho` | Spectral radius | R₀ proxy |

`⊕` is entrywise addition in a fixed partition-, population- and unit-semantics fibre:
the commutative-monoid laws are proved, but addition and categorical coproduct are
distinct constructions and no coproduct universal property is claimed.

`⊗` is Kronecker stratification of a contact matrix, not the slice-category pullback
used to stratify an epidemic model. Model stratification belongs to AlgebraicPetri's
`typed_product`; constructing a contact matrix over a product partition is a separate
operation.

`↓` is not the ordinary left Kan extension of the `MeanContacts` matrix. In the
straightforward discrete, vector-valued interpretation, a left Kan extension along the
partition map — on both contact axes — first produces the coproduct/direct sum over
each fibre. Obtaining a scalar contact total then requires the additive fold, after
which the shipped operation additionally divides by the pushed-forward participant
population.

At the record level, `migrate_coarsen` simply postcomposes participant- and
contact-group assignments with the partition map. For each group-labelled family
`X → G` this is the dependent sum `Σ_f`, equivalently the left Kan extension along the
discrete map `f`; unlike matrix coarsening it needs neither an additive fold nor
population normalisation. A universal property for the migration of the whole
`ContactSurveyACSet` has not been established.

`⊠` coarsens back to the base matrix by design and is regression-tested to, but the
Julia construction validates its mixing distribution and product population only
approximately, so the accepted domain is not literally an exact section; the table
says "assumption-driven" rather than "section of product coarsening" until that is
closed.

## Categorical Framework

### Objects
A `ContactMatrix` bundles:
- An n×n real matrix, read in whichever representation its unit-semantics tag names
- An `AbstractPartition` (age, sex, region, or a product such as age × sex)
- A population vector, required at construction and used by every population-weighted
  operation (symmetrisation, coarsening, refinement, transport, unit conversion)
- Unit semantics (mean contacts / counts / per-capita rates)

### Morphisms
- **Coarsening** (via `PartitionMap`; `AgeMap` is the age-specific alias): partition maps that push contact structure forward along their fibres. Surjectivity is not required — a coarse group hit by no fine group is empty and carries no contacts
- **Refinement** (via `RefinementPrior`): splits each coarse group's contacts between
  its fine subgroups in proportion to a population prior. A section of coarsening —
  coarsening the result recovers the original matrix — but not an inverse, since
  coarsening discards the within-group structure the prior then has to supply
- **Activity refinement** (via `ActivityRefinement`): assumption-driven lift to a product partition such as age × activity
- **Generalized lifts** (via `GeneralizedLift`): assumption-driven lifts from a base partition to a product partition such as age × SES, constructed to coarsen back to the base matrix
- **Constrained partial-data lifts** (via `SourceStratifiedContactMatrix` and `ConstrainedGeneralizedLift`): reconstruct product matrices constrained by participant-side covariates, coarsening, and reciprocity
- **Symmetrisation**: idempotent endomorphism preserving reciprocity
- **Setting composition**: commutative monoid structure (additive)
- **Population transport** (via `transport_population`): changes the population at a
  fixed partition while holding each pair's total contacts fixed. On strictly positive
  populations it satisfies identity, composition and inverse laws, and it commutes
  with coarsening for any partition map. An empty source group whose participant
  column is zero admits a one-way extension of this
- **Demographic reprojection** (via `reproject`): changes the population at a fixed
  partition while establishing reciprocity under the new one, following Arregui et
  al. Unlike population transport it is not functorial in the population and does not
  commute with coarsening — it repairs the matrix rather than carrying it across

### Survey-to-matrix construction
`compute_matrix` (`▷`) normalises the survey's sufficient statistics: the weighted
contact counts `C[i,j]` from participant group `j` into contactee group `i`, and the
weighted participant mass `W_j` per group, giving `M[i,j] = C[i,j]/W_j`.

Building at a fine partition and then coarsening does **not** always agree with
building directly at the coarse partition, because the two normalise by different
quantities: `compute_matrix` divides by participant mass, while `coarsen` weights by
population. Writing `f` for the partition map, the two routes are guaranteed to give
the same **matrix** whenever, within each coarse **participant/source** fibre `J`,

    N_j = λ_J · W_j        for all j ∈ f⁻¹(J)

with `λ_J > 0` on fibres carrying participant mass. The constant cancels in the
derivation, but it must be nonzero; on ContACT's nonnegative population domain this
means strictly positive.

This condition is sufficient for every survey, and is the condition required for
agreement for arbitrary contact-count data. A particular survey may agree without
it — for example, when all fine participant groups in a fibre have the same coarse
contact profile — but the operations themselves do not guarantee it. A fibre with no
participant mass contributes zero source columns on both routes, so the matrix values
agree there for any `N`.

When an external fine population `N` is used, agreement of the matrix values does not
by itself imply agreement of the attached populations. To guarantee agreement of the
full `ContactMatrix` objects, the direct coarse build should also be supplied with the
pushed-forward population `f_*N`.

- With `population=nothing`, the condition holds automatically: `compute_matrix`
  attaches `W` itself as the population, so `λ_J = 1`.
- With an external target population, it becomes a **calibration
  (post-stratification) condition** on the survey weights: weighted participant mass
  must be proportional to the target population within each coarse fibre. The
  proportionality constant may differ between fibres, but the condition must hold in
  every fibre carrying participant mass; supplying `weights=` is not by itself
  sufficient.

When the condition fails, the two routes are generally estimating different quantities
rather than one of them being incorrect, so they need not agree. This behaviour is
pinned by the `"▷ then ↓ versus ▷ at the coarse partition"` testset.

## Formal Proofs (Lean 4)

The `proofs/` directory contains machine-checked proofs of:

| Property | File | Status |
|----------|------|--------|
| Partition maps between contact-matrix objects compose as a category | `ContactCat.lean` | ✅ |
| Category morphisms act functorially on contact data (`coarsenAction`) | `Coarsening.lean` | ✅ |
| Coarsening functoriality (unweighted count model) | `Coarsening.lean` | ✅ |
| Coarsening preserves total contacts (unweighted count model) | `Coarsening.lean` | ✅ |
| **Population-weighted** coarsening functoriality (matches `coarsen`) | `Coarsening.lean` | ✅ |
| **Population-weighted** coarsening preserves total contacts | `Coarsening.lean` | ✅ |
| Refinement is a section of coarsening (empty coarse groups must carry no contacts) | `Refinement.lean` | ✅ |
| That section law with no side condition, where no coarse population is zero | `Refinement.lean` | ✅ |
| Refinement preserves total contacts | `Refinement.lean` | ✅ |
| Refinement inherits reciprocity from the coarse matrix | `Refinement.lean` | ✅ |
| Refinement is **not** a left inverse of coarsening (counterexample) | `Refinement.lean` | ✅ |
| Symmetrisation idempotence (positive populations) | `Symmetrisation.lean` | ✅ |
| Symmetrisation reciprocity (positive populations) | `Symmetrisation.lean` | ✅ |
| Symmetrisation idempotence & reciprocity with empty groups (matches `symmetrise`) | `Symmetrisation.lean` | ✅ |
| Additive composition associativity | `Composition.lean` | ✅ |
| Additive composition commutativity | `Composition.lean` | ✅ |
| Stratification = Kronecker product, with Julia linear-index layout | `Stratification.lean` | ✅ |
| Stratification distributes over composition | `Stratification.lean` | ✅ |
| Symmetrisation–composition commutativity (positive populations) | `Commutativity.lean` | ✅ |
| Symmetrisation–stratification commutativity (iff coupling symmetric) | `Stratification.lean` | ✅ |
| Constrained-lift coarsening recovery, symmetry & structural zeros | `ConstrainedLift.lean` | ✅ |
| Reprojection reciprocity at the target population (incl. empty groups) | `Reprojection.lean` | ✅ |
| Reprojection preserves the grand total of contacts at the target population, summed over all pairs (incl. empty groups) | `Reprojection.lean` | ✅ |
| Reprojection is the identity on matrices already reciprocal at the target | `Reprojection.lean` | ✅ |
| Reprojection is not functorial in the population (counterexample) | `Reprojection.lean` | ✅ |
| Population transport identity & composition, for every matrix (positive populations) | `PopulationTransport.lean` | ✅ |
| Population transport preserves total contacts entrywise (positive target) | `PopulationTransport.lean` | ✅ |
| Population transport preserves reciprocity when the source already has it | `PopulationTransport.lean` | ✅ |
| Population transport is invertible (round trip), for every matrix (positive populations) | `PopulationTransport.lean` | ✅ |
| Population transport commutes with population-weighted coarsening, for any fibre map (nonnegative populations) | `PopulationTransport.lean` | ✅ |
| That commutation genuinely needs a nonzero target population (counterexample) | `PopulationTransport.lean` | ✅ |
| Row/column population rescalings form a groupoid — identity, composition, and inverses where the scalings are nonzero — with the three representations among its objects | `UnitSemantics.lean` | ✅ |
| Conjugating an operation by the representation isomorphisms makes it natural by construction, with no hypothesis on the operation | `UnitSemantics.lean` | ✅ |
| Symmetrisation is natural with respect to the representation isomorphisms | `UnitSemantics.lean` | ✅ |
| Symmetrisation that ignores the input tag is **not** natural (counterexample) | `UnitSemantics.lean` | ✅ |
| Symmetrisation in `ContactCounts` and `PerCapitaRate` is plain averaging | `UnitSemantics.lean` | ✅ |
| Where the effective population is 1 — `ContactCounts` and `PerCapitaRate` — reciprocity is plain matrix symmetry | `UnitSemantics.lean` | ✅ |
| Reciprocal *proportionate* activity lift: build/coarsen equivariance, for any fibre map (nonnegative populations) | `Coarsening.lean` | ✅ |

The equivariance row above is specific to the reciprocal proportionate lift. ContACT's
other two activity kernels, assortative and disassortative, allocate contacts by walking
activity strata in order, and commute with coarsening only when each coarse group is a
contiguous run of fine strata. Merging non-adjacent strata changes the result by O(1) in
mean-contacts units — a documented limitation pinned by regression test rather than a
proof. Ordinary use, merging adjacent strata of an ordered activity partition, satisfies
the condition automatically.

The weighted-coarsening and empty-group theorems model the exact arithmetic the
Julia code runs (population-weighted averaging and the zero-population branch),
not only the simplified count-space model; the count-space theorems are retained
because the weighted results are built on top of them.

The population-transport theorems are equalities in ℝ. In floating point,
`transport_population`'s identity law is bitwise exact — a column whose population
is unchanged is copied rather than rescaled, since `M·N/N` is not the identity in
Float64 — while composition, invertibility and total-contact preservation agree to
roundoff, about one ulp. Its Julia domain is also narrower than the theorems':
a zero source population is admitted only when that group's **participant column** is
zero. Its row may be nonzero, and which axis is guarded matters — the corresponding
guards elsewhere in the package check different axes.

The `ContactCat.lean` row covers less than its name suggests: its morphism carries a
map between index sets and a surjectivity proof, so the category laws are about those
maps rather than about the contact data. The theorems that do relate a coarsened matrix
to the one it came from are the coarsening rows above.

Build proofs with:
```bash
cd proofs && lake build
```

## Installation

ContACT is not yet on the General registry, so install it from the repository:

```julia
using Pkg
Pkg.add(url="https://github.com/epirecipes/ContACT.jl")
```

## Dependencies

- [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl) — categorical algebra (FinFunctions for partition maps)
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) — survey data handling
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) — matrix operations
- [Random](https://docs.julialang.org/en/v1/stdlib/Random/) — sampling and MCMC over constrained lifts

## Related Packages

- [socialmixr](https://github.com/epiforecasts/socialmixr) (R) — the original contact matrix toolkit
- [CategoricalPopulationDynamics.jl](https://github.com/ecorecipes/CategoricalPopulationDynamics.jl) — categorical population models (shares coarsening/stratification machinery)
