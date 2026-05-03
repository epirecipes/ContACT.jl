# Introduction to ContACT.jl
Simon Frost
2026-05-03

- [Overview](#overview)
- [Setup](#setup)
- [Loading POLYMOD Data](#loading-polymod-data)
- [Preparing the Survey](#preparing-the-survey)
- [Computing a Contact Matrix](#computing-a-contact-matrix)
- [Symmetrisation (Idempotent
  Endomorphism)](#symmetrisation-idempotent-endomorphism)
- [Coarsening (Left Kan Extension)](#coarsening-left-kan-extension)
  - [Functoriality](#functoriality)
- [What’s Next?](#whats-next)

## Overview

ContACT.jl applies **category theory** to the construction and
manipulation of age-structured contact matrices from social mixing
surveys. This vignette introduces the core workflow using the POLYMOD
study (Mossong et al., 2008) — the gold standard for parameterising
infectious disease models.

The key categorical insight: a contact matrix is not merely a numerical
array but an **object** in a category where morphisms are
structure-preserving transformations (coarsening, symmetrisation,
composition).

## Setup

``` julia
using ContACT
using CSV
using DataFrames
using LinearAlgebra
```

## Loading POLYMOD Data

We use a subset of the POLYMOD UK study: 1012 participants and their
11876 reported contacts.

``` julia
data_dir = joinpath(@__DIR__, "..", "..", "data")
participants = CSV.read(joinpath(data_dir, "polymod_uk_participants.csv"), DataFrame)
contacts = CSV.read(joinpath(data_dir, "polymod_uk_contacts.csv"), DataFrame)

println("Participants: $(nrow(participants))")
println("Contacts: $(nrow(contacts))")
first(participants, 5)
```

    Participants: 1012
    Contacts: 11876

<div><div style = "float: left;"><span>5×4 DataFrame</span></div><div style = "clear: both;"></div></div><div class = "data-frame" style = "overflow-x: scroll;">

| Row | part_id | part_age_exact | country        | dayofweek |
|----:|--------:|---------------:|:---------------|:----------|
|     |   Int64 |          Int64 | String15       | String3   |
|   1 |    4517 |              5 | United Kingdom | 4         |
|   2 |    4518 |              5 | United Kingdom | 4         |
|   3 |    4519 |              6 | United Kingdom | 4         |
|   4 |    4520 |              3 | United Kingdom | 6         |
|   5 |    4521 |              2 | United Kingdom | 6         |

</div>

## Preparing the Survey

ContACT.jl expects a `ContactSurvey` with columns `:part_id`,
`:part_age` (participants) and `:part_id`, `:cnt_age` (contacts).
POLYMOD stores ages as `:part_age_exact` and either `:cnt_age_exact` or
estimated ranges (as strings with “NA” for missing).

``` julia
# Helper to parse string ages (CSV stores them as strings with "NA")
parse_age(x::AbstractString) = x == "NA" ? missing : parse(Float64, x)
parse_age(x::Real) = Float64(x)
parse_age(::Missing) = missing

# Rename participant age
rename!(participants, :part_age_exact => :part_age)
participants.part_age = Float64.(participants.part_age)

# Resolve contact ages: use exact if available, otherwise midpoint of range
cnt_exact = parse_age.(contacts.cnt_age_exact)
cnt_min = parse_age.(contacts.cnt_age_est_min)
cnt_max = parse_age.(contacts.cnt_age_est_max)
contacts.cnt_age = coalesce.(cnt_exact, (cnt_min .+ cnt_max) ./ 2)
select!(contacts, [:part_id, :cnt_age])
dropmissing!(contacts, :cnt_age)

survey = ContactSurvey(participants, contacts)
println("Survey: $(nrow(survey.participants)) participants, $(nrow(survey.contacts)) contacts")
```

    Survey: 1012 participants, 11873 contacts

## Computing a Contact Matrix

The central operation is the **restricted functor** from surveys to
contact matrices. We fix an age partition — this makes the operation
deterministic and functorial:

``` julia
# Standard 5-year age bands up to 75+
partition = AgePartition(
    [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
)

# UK 2005 population (thousands) — approximation from ONS
uk_pop = [3433.0, 3554.0, 3824.0, 3956.0, 3757.0, 3520.0, 4009.0,
           4405.0, 4548.0, 4187.0, 3883.0, 3589.0, 3117.0, 2765.0,
           2425.0, 3752.0]

# The ▷ operator applies the survey→matrix functor (type \triangleright<TAB>)
cm = survey ▷ partition
println("Contact matrix: $(n_groups(cm)) age groups")
println("Spectral radius ρ(M) ∝ R₀: $(round(ρ(cm); digits=2))")
```

    Contact matrix: 16 age groups
    Spectral radius ρ(M) ∝ R₀: 12.24

The result is a `ContactMatrix` — not just a bare array, but a
categorical object bundling the matrix with its age partition,
population, and semantics:

``` julia
println("Type: $(typeof(cm))")
println("Semantics: $(cm.semantics)")
println("Age groups: $(age_labels(cm))")
```

    Type: ContACT.ContactMatrix{Float64, ContACT.MeanContacts}
    Semantics: ContACT.MeanContacts()
    Age groups: ["[0,5)", "[5,10)", "[10,15)", "[15,20)", "[20,25)", "[25,30)", "[30,35)", "[35,40)", "[40,45)", "[45,50)", "[50,55)", "[55,60)", "[60,65)", "[65,70)", "[70,75)", "75+"]

## Symmetrisation (Idempotent Endomorphism)

Raw survey matrices are asymmetric because group sizes differ. The
reciprocity constraint ensures total contacts from group $i$ to $j$
equal total contacts from $j$ to $i$:

$$M_{ij} \cdot N_j = M_{ji} \cdot N_i$$

Symmetrisation is an **idempotent endomorphism** in the contact matrix
category. Use `↔` (`\leftrightarrow<TAB>`) for the reciprocity
projection:

``` julia
cm_sym = ↔(cm)

# Verify reciprocity
M = matrix(cm_sym)
N = population(cm_sym)
max_violation = maximum(abs(M[i,j] * N[j] - M[j,i] * N[i])
    for i in 1:n_groups(cm_sym), j in 1:n_groups(cm_sym))
println("Max reciprocity violation: $(max_violation)")

# Idempotence: applying twice gives the same result
cm_sym2 = ↔(cm_sym)
println("Idempotent: $(matrix(cm_sym2) ≈ matrix(cm_sym))")
```

    Max reciprocity violation: 1.4210854715202004e-14
    Idempotent: true

## Coarsening (Left Kan Extension)

A key morphism is **coarsening**: reducing the number of age groups.
Categorically, this is a left Kan extension along a surjective age-group
map $f: G_{\mathrm{fine}} \to G_{\mathrm{coarse}}$.

The coarse partition must have limits that are a **subset** of the fine
partition’s limits (ensuring a well-defined surjection between groups).

``` julia
# Coarsen to broader groups (limits must be a subset of the fine partition)
coarse = AgePartition([0, 15, 45, 65])
cm_coarse = cm ↓ coarse

println("Coarsened: $(n_groups(cm_coarse)) groups")
println("Age groups: $(age_labels(cm_coarse))")
println("Matrix:")
display(round.(matrix(cm_coarse); digits=2))
```

    Coarsened: 4 groups
    Age groups: ["[0,15)", "[15,45)", "[45,65)", "65+"]
    Matrix:

    4×4 Matrix{Float64}:
     6.85  2.2   0.98  0.89
     4.39  7.6   4.77  3.25
     1.33  2.27  3.18  2.42
     0.27  0.46  0.85  1.72

### Functoriality

The defining property of a functor: coarsening in two steps gives the
same result as coarsening in one step. Using `∘` for map composition:

``` julia
# Two-step: 16 → 4 → 2
mid = AgePartition([0, 15, 45, 65])
final_part = AgePartition([0, 45])

# Compose the maps: fine → mid → coarse
f = AgeMap(partition, mid)
g = AgeMap(mid, final_part)
h = g ∘ f   # composed map (type \circ<TAB>)

# Functoriality: (cm ↓ g) ∘ (cm ↓ f) == cm ↓ (g ∘ f)
via_mid = (cm ↓ f) ↓ g
direct = cm ↓ h

println("Functoriality: cm ↓ (g ∘ f) == (cm ↓ f) ↓ g")
println("  Holds: $(matrix(via_mid) ≈ matrix(direct))")
println("\nVia mid:")
display(round.(matrix(via_mid); digits=3))
println("\nDirect:")
display(round.(matrix(direct); digits=3))
```

    Functoriality: cm ↓ (g ∘ f) == (cm ↓ f) ↓ g
      Holds: true

    Via mid:

    Direct:

    2×2 Matrix{Float64}:
     10.406  5.436
      2.262  4.054

    2×2 Matrix{Float64}:
     10.406  5.436
      2.262  4.054

## What’s Next?

- **Vignette 2**: Composition and stratification — combining settings
  and spatial structure
- **Vignette 3**: The categorical framework — ACSets, schemas, and
  functorial data migration
