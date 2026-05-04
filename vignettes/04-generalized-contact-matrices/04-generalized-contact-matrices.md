# Generalized Contact Matrices
Simon Frost
2026-05-03

- [Overview](#overview)
- [Setup](#setup)
- [A Base Age Matrix](#a-base-age-matrix)
- [Add an SES Dimension](#add-an-ses-dimension)
- [Random/Product Mixing Lift](#randomproduct-mixing-lift)
- [Assortative SES Mixing](#assortative-ses-mixing)
- [Reciprocity and Total Contacts](#reciprocity-and-total-contacts)
- [R₀ Utilities](#r₀-utilities)
- [Custom Kernels](#custom-kernels)
- [Interpretation](#interpretation)

## Overview

Manna et al. describe **generalized contact matrices** whose rows and
columns are indexed not only by age, but by product groups such as age ×
wealth, age × education, or age × ethnicity. In ContACT.jl this is not a
separate matrix type: it is a `ContactMatrix` over a `ProductPartition`.

If a survey records both participant and contact-side SES variables, the
survey-to-matrix functor can estimate the product matrix directly:

``` julia
G = survey ▷ (age × wealth)
```

When the contact-side SES variable is not observed, the generalized
matrix is an assumption-driven **lift** from a base partition `P` to a
product partition `P × Q`. The lift is a section of the projection
`P × Q → P`: coarsening back to `P` recovers the original matrix.

## Setup

``` julia
using ContACT
using LinearAlgebra

import ContACT: ⊠, ↓, ↔, ρ, ×
```

## A Base Age Matrix

Start with an age-structured contact matrix. The entries are mean
contacts made with the row group by participants in the column group. We
symmetrise so total contacts are reciprocal before lifting.

``` julia
age = AgePartition([0, 18, 65]; labels=["child", "adult", "older"])
age_pop = [1200.0, 3600.0, 1400.0]

M_age = [3.5 0.8 0.2;
         1.2 4.0 0.7;
         0.4 1.4 2.0]

cm_age = ↔(ContactMatrix(M_age, age, age_pop))
round.(matrix(cm_age); digits=2)
```

    3×3 Matrix{Float64}:
     3.5   0.6   0.27
     1.8   4.0   2.15
     0.32  0.84  2.0

## Add an SES Dimension

The second dimension is an ordinary categorical partition. It does not
need to be age-specific; the same machinery works for wealth, education,
ethnicity, occupation, risk class, or any finite grouping.

``` julia
wealth = CategoricalPartition(:wealth;
    levels=["low", "middle", "high"],
    labels=["low SES", "middle SES", "high SES"],
)

wealth_distribution = [0.35, 0.45, 0.20]
group_labels(age × wealth)
```

    9-element Vector{String}:
     "child:low SES"
     "child:middle SES"
     "child:high SES"
     "adult:low SES"
     "adult:middle SES"
     "adult:high SES"
     "older:low SES"
     "older:middle SES"
     "older:high SES"

## Random/Product Mixing Lift

Random mixing splits each age-age total-contact block in proportion to
the product of the SES distribution. The operator `⊠` (`\boxtimes<TAB>`)
denotes the product lift.

``` julia
random_spec = GeneralizedLift(wealth; distribution=wealth_distribution)
G_random = cm_age ⊠ random_spec

println("Groups: $(n_groups(G_random))")
println("Coarsens back to age matrix: $(matrix(G_random ↓ age) ≈ matrix(cm_age))")
println("ρ(base)   = $(round(ρ(cm_age); digits=3))")
println("ρ(random) = $(round(ρ(G_random); digits=3))")
```

    Groups: 9
    Coarsens back to age matrix: true
    ρ(base)   = 5.301
    ρ(random) = 5.301

Categorically, the product projection is a partition morphism:

``` julia
π_age = PartitionMap(G_random.partition, age)
matrix(G_random ↓ π_age) ≈ matrix(cm_age)
```

    true

## Assortative SES Mixing

Manna et al. also consider an assortative construction in which the
added dimension has its own activity shares and within-group contact
fractions. For a three-level SES dimension, ContACT.jl provides
`AssortativeDimensionMixing` with the same parameterization:

``` julia
ses_activity = [0.20, 0.40, 0.40]
ses_assortativity = [0.60, 0.50, 0.65]
offdiag_split = [0.60, 0.60, 0.50]

assortative = AssortativeDimensionMixing(
    ses_activity,
    ses_assortativity;
    offdiag_split=offdiag_split,
)

assortative_spec = GeneralizedLift(
    wealth;
    distribution=wealth_distribution,
    mixing=assortative,
)

G_assortative = cm_age ⊠ assortative_spec
println("Coarsens back to age matrix: $(matrix(G_assortative ↓ age) ≈ matrix(cm_age))")
println("ρ(assortative) = $(round(ρ(G_assortative); digits=3))")
```

    Coarsens back to age matrix: true
    ρ(assortative) = 7.815

The age-only matrix is recovered by coarsening, but the full matrix
contains additional within-age SES structure:

``` julia
round.(matrix(G_assortative)[1:3, 1:3]; digits=2)
```

    3×3 Matrix{Float64}:
     1.2  0.54  0.17
     0.7  1.56  2.27
     0.1  1.01  4.55

We can also marginalize to the SES-only matrix:

``` julia
G_wealth = G_assortative ↓ wealth
round.(matrix(G_wealth); digits=2)
```

    3×3 Matrix{Float64}:
     1.8   0.79  0.63
     1.02  2.33  3.09
     0.36  1.38  6.81

## Reciprocity and Total Contacts

Because the base matrix is reciprocal and the lift uses symmetric
diagonal blocks with transposed off-diagonal blocks, reciprocity is
preserved in total-contact space.

``` julia
totals = matrix(G_assortative) * Diagonal(population(G_assortative))
maximum(abs.(totals - transpose(totals)))
```

    5.684341886080802e-14

## R₀ Utilities

The paper compares epidemic thresholds under age-only, SES-only, and
generalized matrices. ContACT.jl provides a next-generation matrix
convention and `R₀` helpers for this calculation.

``` julia
γ = 0.25
β = calibrate_transmissibility(cm_age, 2.7; recovery_rate=γ)

r0_values = [
    "age only" => R₀(cm_age; transmissibility=β, recovery_rate=γ),
    "age × wealth, random" => R₀(G_random; transmissibility=β, recovery_rate=γ),
    "age × wealth, assortative" => R₀(G_assortative; transmissibility=β, recovery_rate=γ),
    "wealth only marginal" => R₀(G_wealth; transmissibility=β, recovery_rate=γ),
]

[(label, round(value; digits=3)) for (label, value) in r0_values]
```

    4-element Vector{Tuple{String, Float64}}:
     ("age only", 2.7)
     ("age × wealth, random", 2.7)
     ("age × wealth, assortative", 3.981)
     ("wealth only marginal", 3.932)

## Custom Kernels

For dimensions or assumptions not covered by the built-in random and
three-level assortative constructors, use `BlockMixing` directly. The
diagonal block is used within a base group; the off-diagonal block is
used for one direction and transposed for the reverse direction.

``` julia
block = wealth_distribution * transpose(wealth_distribution)
custom_spec = GeneralizedLift(
    wealth;
    distribution=wealth_distribution,
    mixing=BlockMixing(block),
)

matrix((cm_age ⊠ custom_spec) ↓ age) ≈ matrix(cm_age)
```

    true

## Interpretation

The generalized lift is not an inverse to coarsening. It is a chosen
section:

``` math
P \times Q \xrightarrow{\;\pi_P\;} P
```

Coarsening along $\pi_P$ is canonical and information-losing. Lifting
from `P` to `P × Q` requires extra data: a population distribution and a
mixing kernel. This is exactly the categorical distinction between a
functorial pushforward and an assumption-driven reconstruction.
