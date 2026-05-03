# Composition and Stratification
Simon Frost
2026-05-03

- [Overview](#overview)
- [Setup](#setup)
- [Setting-Specific Matrices](#setting-specific-matrices)
- [Additive Composition with ⊕](#additive-composition-with-)
- [Intervention Modelling](#intervention-modelling)
- [Spatial Stratification with ⊗](#spatial-stratification-with-)
- [Composing Operations](#composing-operations)
- [What’s Next?](#whats-next)

## Overview

Contact surveys record the **setting** where contacts occur: home, work,
school, transport, leisure. The total contact matrix decomposes as a sum
of setting-specific matrices. This is not merely addition — it’s a
**commutative monoid** structure in the category of contact matrices
with formal guarantees:

- **Associativity**: $(A \oplus B) \oplus C = A \oplus (B \oplus C)$
- **Commutativity**: $A \oplus B = B \oplus A$
- **Identity**: $A \oplus \mathbf{0} = A$

Spatial stratification uses the **Kronecker product** ($\otimes$) with a
coupling matrix to create multi-region contact patterns.

## Setup

``` julia
using ContACT
using CSV
using DataFrames
using LinearAlgebra
```

## Setting-Specific Matrices

We construct synthetic setting-specific matrices based on known patterns
from POLYMOD for a 4-group age structure:

``` julia
partition = AgePartition([0, 5, 18, 65])
uk_pop = [6987.0, 11537.0, 35854.0, 9492.0]  # aggregated to 4 groups (thousands)

# Home: strong assortative mixing, parent-child contacts
M_home = [1.8 0.3 0.4 0.1;
          0.3 0.6 0.2 0.1;
          0.4 0.2 1.0 0.3;
          0.1 0.1 0.3 0.8]

# Work: concentrated in working-age adults
M_work = [0.0 0.0 0.0 0.0;
          0.0 0.1 1.5 0.2;
          0.0 1.5 2.5 0.4;
          0.0 0.2 0.4 0.3]

# School: strong child-child contacts
M_school = [3.5 2.0 0.1 0.0;
            2.0 5.0 0.3 0.0;
            0.1 0.3 0.2 0.0;
            0.0 0.0 0.0 0.0]

# Other: weak mixing across all groups
M_other = [0.5 0.3 0.4 0.2;
           0.3 0.5 0.8 0.3;
           0.4 0.8 1.2 0.5;
           0.2 0.3 0.5 0.4]

cm_home   = ContactMatrix(M_home, partition, uk_pop)
cm_work   = ContactMatrix(M_work, partition, uk_pop)
cm_school = ContactMatrix(M_school, partition, uk_pop)
cm_other  = ContactMatrix(M_other, partition, uk_pop)
```

    ContactMatrix{Float64, MeanContacts} (4×4 age groups)

## Additive Composition with ⊕

The total contact matrix is the **monoidal sum** of all settings:

``` julia
# Using the ⊕ operator (type \oplus<TAB> in Julia REPL)
cm_total = cm_home ⊕ cm_work ⊕ cm_school ⊕ cm_other

println("Total contact matrix (4 age groups):")
display(round.(matrix(cm_total); digits=2))
println("\nSpectral radius ρ(M): $(round(ρ(cm_total); digits=2))")
```

    Total contact matrix (4 age groups):

    Spectral radius ρ(M): 10.21

    4×4 Matrix{Float64}:
     5.8  2.6  0.9  0.3
     2.6  6.2  2.8  0.6
     0.9  2.8  4.9  1.2
     0.3  0.6  1.2  1.5

Verify the monoid axioms:

``` julia
# Associativity
lhs = (cm_home ⊕ cm_work) ⊕ cm_school
rhs = cm_home ⊕ (cm_work ⊕ cm_school)
println("Associativity: $(matrix(lhs) ≈ matrix(rhs))")

# Commutativity
println("Commutativity: $(matrix(cm_home ⊕ cm_work) ≈ matrix(cm_work ⊕ cm_home))")

# Identity (zero matrix)
zero_cm = ContactMatrix(zeros(4, 4), partition, uk_pop)
println("Identity: $(matrix(cm_home ⊕ zero_cm) ≈ matrix(cm_home))")
```

    Associativity: true
    Commutativity: true
    Identity: true

## Intervention Modelling

The compositional structure makes it natural to model interventions that
selectively modify settings:

``` julia
# School closure: remove school contacts entirely
cm_closure = cm_home ⊕ cm_work ⊕ cm_other

# Work-from-home: reduce work contacts by 70%
M_work_reduced = 0.3 .* M_work
cm_wfh = ContactMatrix(M_work_reduced, partition, uk_pop)
cm_lockdown = cm_home ⊕ cm_wfh ⊕ cm_other

# Full lockdown: only home contacts
cm_full_lockdown = cm_home

# Compare spectral radii (proportional to R₀)
println("Scenario analysis — ρ(M) as R₀ proxy:")
println("  Normal:          $(round(ρ(cm_total); digits=2))")
println("  School closure:  $(round(ρ(cm_closure); digits=2))")
println("  Lockdown (WFH):  $(round(ρ(cm_lockdown); digits=2))")
println("  Full lockdown:   $(round(ρ(cm_full_lockdown); digits=2))")
```

    Scenario analysis — ρ(M) as R₀ proxy:
      Normal:          10.21
      School closure:  6.62
      Lockdown (WFH):  4.59
      Full lockdown:   2.09

## Spatial Stratification with ⊗

When modelling epidemics across regions, we combine **local contact
patterns** with **inter-regional coupling**. This is the Kronecker
product:

$$M_{\text{spatial}} = C \otimes M_{\text{local}}$$

where $C$ is the coupling matrix encoding between-region mixing.

``` julia
# Three regions with mostly-local mixing
coupling = [0.8 0.15 0.05;
            0.15 0.7 0.15;
            0.05 0.15 0.8]

# Stratify the total contact matrix
cm_spatial = cm_total ⊗ coupling

println("Stratified matrix: $(n_groups(cm_spatial)) groups (3 regions × 4 ages)")
println("ρ(M_spatial): $(round(ρ(cm_spatial); digits=2))")
```

    Stratified matrix: 12 groups (3 regions × 4 ages)
    ρ(M_spatial): 10.21

The resulting block structure:

``` julia
M_spatial = matrix(cm_spatial)
println("Block structure (region 1 ↔ region 1, top-left 4×4):")
display(round.(M_spatial[1:4, 1:4]; digits=2))
println("\nBlock structure (region 1 ↔ region 2, off-diagonal 4×4):")
display(round.(M_spatial[1:4, 5:8]; digits=2))
```

    Block structure (region 1 ↔ region 1, top-left 4×4):

    Block structure (region 1 ↔ region 2, off-diagonal 4×4):

    4×4 Matrix{Float64}:
     4.64  2.08  0.72  0.24
     2.08  4.96  2.24  0.48
     0.72  2.24  3.92  0.96
     0.24  0.48  0.96  1.2

    4×4 Matrix{Float64}:
     0.87  0.39  0.14  0.05
     0.39  0.93  0.42  0.09
     0.14  0.42  0.74  0.18
     0.05  0.09  0.18  0.22

## Composing Operations

The categorical framework ensures operations compose cleanly. Using
ContACT.jl’s Unicode operators, the full pipeline reads algebraically:

``` julia
# Start with fine-grained POLYMOD data
data_dir = joinpath(@__DIR__, "..", "..", "data")
participants = CSV.read(joinpath(data_dir, "polymod_uk_participants.csv"), DataFrame)
contacts_df = CSV.read(joinpath(data_dir, "polymod_uk_contacts.csv"), DataFrame)

parse_age(x::AbstractString) = x == "NA" ? missing : parse(Float64, x)
parse_age(x::Real) = Float64(x)
parse_age(::Missing) = missing

rename!(participants, :part_age_exact => :part_age)
participants.part_age = Float64.(participants.part_age)

cnt_exact = parse_age.(contacts_df.cnt_age_exact)
cnt_min = parse_age.(contacts_df.cnt_age_est_min)
cnt_max = parse_age.(contacts_df.cnt_age_est_max)
contacts_df.cnt_age = coalesce.(cnt_exact, (cnt_min .+ cnt_max) ./ 2)
select!(contacts_df, [:part_id, :cnt_age])
dropmissing!(contacts_df, :cnt_age)

survey = ContactSurvey(participants, contacts_df)

# ─── The algebraic pipeline ─────────────────────────────────────────────
# Each line uses a Unicode operator:
#   ▷  functor application (survey → matrix)
#   ↓  coarsening (left Kan extension)
#   ⊗  stratification (Kronecker product)
#   ρ  spectral radius

fine_part = AgePartition([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75])
cm_fine = survey ▷ fine_part                    # F(survey)
cm_4grp = cm_fine ↓ AgePartition([0, 15, 45, 65])   # left Kan extension
cm_sym  = ↔(cm_4grp)                            # reciprocity projection
cm_3reg = cm_sym ⊗ coupling                     # Kronecker stratification

println("Pipeline: Survey ▷ Partition ↓ Coarse ↔ ⊗ Coupling")
println("  $(n_groups(cm_fine)) → $(n_groups(cm_4grp)) → $(n_groups(cm_3reg)) groups")
println("  ρ(M_fine)   = $(round(ρ(cm_fine); digits=2))")
println("  ρ(M_coarse) = $(round(ρ(cm_4grp); digits=2))")
println("  ρ(M_sym)    = $(round(ρ(cm_sym); digits=2))")
println("  ρ(M_3reg)   = $(round(ρ(cm_3reg); digits=2))")
```

    Pipeline: Survey ▷ Partition ↓ Coarse ↔ ⊗ Coupling
      16 → 4 → 12 groups
      ρ(M_fine)   = 12.24
      ρ(M_coarse) = 11.92
      ρ(M_sym)    = 11.98
      ρ(M_3reg)   = 11.98

## What’s Next?

- **Vignette 3**: The categorical framework — ACSets, schemas,
  functorial data migration, and undirected wiring diagrams for
  declarative composition
