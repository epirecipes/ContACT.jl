# ContACT.jl

**Category theory applied to contact matrices** — a Julia package for building, composing, and manipulating structured contact matrices using applied category theory.

## Overview

ContACT.jl provides a categorical framework for:
- Computing contact matrices from individual-level survey data over arbitrary finite partitions (the **functor** from surveys to matrices)
- Composing setting-specific matrices (home, work, school) via **additive composition** (⊕)
- Changing partition resolution via **coarsening** (left Kan extension) and refinement
- Refining reciprocal matrices by respondent social activity, following Britton-Ball-style high/low or quantile activity strata (⤊)
- Lifting age matrices to generalized age × SES/contact matrices with explicit mixing assumptions (⊠)
- Reconstructing generalized matrices from partial participant-side covariates with source-stratified intermediate matrices
- Adding spatial structure via **stratification** (Kronecker product) (⊗)
- Enforcing reciprocity via **symmetrisation** (↔)

All operations come with formal guarantees verified in Lean 4.

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

| Symbol | LaTeX | Operation | Categorical Meaning |
|--------|-------|-----------|-------------------|
| `⊕` | `\oplus` | Additive composition | Coproduct in setting category |
| `⊗` | `\otimes` | Stratification | Pullback in slice category |
| `↓` | `\downarrow` | Coarsening | Left Kan extension |
| `↑` | `\uparrow` | Refinement with prior | Parameterised disaggregation |
| `⤊` | `\Uuparrow` | Activity refinement | Hidden-stratum lift |
| `⊠` | `\boxtimes` | Generalized product lift | Section of product coarsening |
| `▷` | `\triangleright` | Survey-to-matrix functor | Functor application |
| `∘` | `\circ` | PartitionMap composition | Morphism composition |
| `↔` | `\leftrightarrow` | Symmetrisation | Reciprocity projection |
| `ρ` | `\rho` | Spectral radius | R₀ proxy |

## Categorical Framework

### Objects
A `ContactMatrix` bundles:
- An n×n real matrix of mean contacts
- An `AbstractPartition` (age, sex, region, or a product such as age × sex)
- A population vector (required for symmetrisation)
- Unit semantics (mean contacts / counts / per-capita rates)

### Morphisms
- **Coarsening** (via `PartitionMap`; `AgeMap` is the age-specific alias): surjective partition maps that push forward contact structure
- **Activity refinement** (via `ActivityRefinement`): assumption-driven lift to a product partition such as age × activity
- **Generalized lifts** (via `GeneralizedLift`): assumption-driven sections from a base partition to a product partition such as age × SES
- **Constrained partial-data lifts** (via `SourceStratifiedContactMatrix` and `ConstrainedGeneralizedLift`): reconstruct product matrices constrained by participant-side covariates, coarsening, and reciprocity
- **Symmetrisation**: idempotent endomorphism preserving reciprocity
- **Setting composition**: commutative monoid structure (additive)

### Functor
`compute_matrix` is a functor from the subcategory of surveys (with fixed partition and weighting) to the category of contact matrices.

## Formal Proofs (Lean 4)

The `proofs/` directory contains machine-checked proofs of:

| Property | File | Status |
|----------|------|--------|
| Contact matrices form a category | `ContactCat.lean` | ✅ |
| Category morphisms act functorially on contact data (`coarsenAction`) | `Coarsening.lean` | ✅ |
| Coarsening functoriality (unweighted count model) | `Coarsening.lean` | ✅ |
| Coarsening preserves total contacts (unweighted count model) | `Coarsening.lean` | ✅ |
| **Population-weighted** coarsening functoriality (matches `coarsen`) | `Coarsening.lean` | ✅ |
| **Population-weighted** coarsening preserves total contacts | `Coarsening.lean` | ✅ |
| Symmetrisation idempotence (positive populations) | `Symmetrisation.lean` | ✅ |
| Symmetrisation reciprocity (positive populations) | `Symmetrisation.lean` | ✅ |
| Symmetrisation idempotence & reciprocity with empty groups (matches `symmetrise`) | `Symmetrisation.lean` | ✅ |
| Additive composition associativity | `Composition.lean` | ✅ |
| Additive composition commutativity | `Composition.lean` | ✅ |
| Stratification = Kronecker product, with Julia linear-index layout | `Stratification.lean` | ✅ |
| Stratification distributes over composition | `Stratification.lean` | ✅ |
| Symmetrisation–composition commutativity | `Commutativity.lean` | ✅ |
| Symmetrisation–stratification commutativity (iff coupling symmetric) | `Stratification.lean` | ✅ |
| Constrained-lift coarsening recovery, symmetry & structural zeros | `ConstrainedLift.lean` | ✅ |

The weighted-coarsening and empty-group theorems model the exact arithmetic the
Julia code runs (population-weighted averaging and the zero-population branch),
not only the simplified count-space model; the count-space theorems are retained
because the weighted results are built on top of them.

Build proofs with:
```bash
cd proofs && lake build
```

## Installation

```julia
using Pkg
Pkg.add("ContACT")
```

For the development version:
```julia
Pkg.add(url="https://github.com/epirecipes/ContACT.jl")
```

## Dependencies

- [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl) — categorical algebra (FinFunctions for partition maps)
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl) — survey data handling
- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) — matrix operations

## Related Packages

- [socialmixr](https://github.com/epiforecasts/socialmixr) (R) — the original contact matrix toolkit
- [CategoricalPopulationDynamics.jl](https://github.com/ecorecipes/CategoricalPopulationDynamics.jl) — categorical population models (shares coarsening/stratification machinery)
