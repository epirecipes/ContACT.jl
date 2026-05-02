# The Categorical Framework
Simon Frost
2026-05-02

- [Overview](#overview)
- [Setup](#setup)
- [ACSet Schemas](#acset-schemas)
  - [ContactSurveySchema](#contactsurveyschema)
  - [ContactMatrixSchema](#contactmatrixschema)
- [Functorial Data Migration](#functorial-data-migration)
  - [Functoriality of Migration](#functoriality-of-migration)
- [Undirected Wiring Diagrams for
  Composition](#undirected-wiring-diagrams-for-composition)
  - [Intervention Scenarios via Diagram
    Modification](#intervention-scenarios-via-diagram-modification)
- [Putting It All Together](#putting-it-all-together)
- [Summary](#summary)

## Overview

This vignette dives into the **deeper categorical machinery** of
ContACT.jl, powered by
[Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl). While the
previous vignettes used the simplified API (plain structs + operators),
the ACSet layer provides:

1.  **Schema-enforced structure**: Type-safe relationships between
    participants, contacts, and age groups, preventing structural errors
    at construction.
2.  **Functorial data migration** ($\Sigma_f$, $\Delta_f$): Principled
    coarsening and restriction of age-structured data along morphisms.
3.  **Undirected wiring diagrams (UWDs)**: Declarative specification of
    how setting-specific matrices compose via shared boundaries.

## Setup

``` julia
using ContACT
using CSV
using DataFrames
using LinearAlgebra
using Catlab
using Catlab.CategoricalAlgebra
import ContACT: ⊕, ⊗, ↓
```

## ACSet Schemas

An **Attributed C-Set** (ACSet) is a functor from a schema category to
**Set**. ContACT.jl defines two schemas:

### ContactSurveySchema

The survey schema encodes the relational structure:

- **Objects**: P (Participants), C (Contacts), G (Age Groups)
- **Morphisms**: `part_group: P → G`, `cnt_group: C → G`,
  `reporter: C → P`

This makes the foreign-key relationships explicit and type-checked.

``` julia
# Load POLYMOD data
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
```

    ContactSurvey(1012×4 DataFrame
      Row │ part_id  part_age  country         dayofweek 
          │ Int64    Float64   String15        String3   
    ──────┼──────────────────────────────────────────────
        1 │    4517       5.0  United Kingdom  4
        2 │    4518       5.0  United Kingdom  4
        3 │    4519       6.0  United Kingdom  4
        4 │    4520       3.0  United Kingdom  6
        5 │    4521       2.0  United Kingdom  6
        6 │    4522       4.0  United Kingdom  6
        7 │    4523       6.0  United Kingdom  6
        8 │    4524       7.0  United Kingdom  6
      ⋮   │    ⋮        ⋮            ⋮             ⋮
     1006 │    5516      63.0  United Kingdom  5
     1007 │    5517      59.0  United Kingdom  5
     1008 │    5518      50.0  United Kingdom  5
     1009 │    5519      57.0  United Kingdom  5
     1010 │    5520      52.0  United Kingdom  5
     1011 │    5521      34.0  United Kingdom  5
     1012 │    5522      39.0  United Kingdom  5
                                         997 rows omitted, 11873×2 DataFrame
       Row │ part_id  cnt_age 
           │ Int64    Float64 
    ───────┼──────────────────
         1 │    4517      4.0
         2 │    4517     40.0
         3 │    4517     31.0
         4 │    4517     52.5
         5 │    4517     29.0
         6 │    4517     59.0
         7 │    4517     57.0
         8 │    4517      4.0
       ⋮   │    ⋮        ⋮
     11867 │    5521      0.0
     11868 │    5521     40.0
     11869 │    5522     15.0
     11870 │    5522     35.0
     11871 │    5522     50.0
     11872 │    5522     35.0
     11873 │    5522     45.0
            11858 rows omitted, Dict{Symbol, Any}())

Convert to an ACSet representation:

``` julia
partition = AgePartition([0, 5, 18, 45, 65])
acs = ContactSurveyACSet(survey, partition)

println("Survey ACSet:")
println("  Age groups (G): $(nparts(acs, :G))")
println("  Participants (P): $(nparts(acs, :P))")
println("  Contacts (C): $(nparts(acs, :C))")
```

    Survey ACSet:
      Age groups (G): 5
      Participants (P): 1012
      Contacts (C): 11873

The ACSet enforces referential integrity: every contact’s `:reporter`
points to a valid participant, and every `:part_group`/`:cnt_group`
points to a valid age group.

``` julia
# Check: all reporters are valid participants
all_valid = all(1 .<= subpart(acs, :reporter) .<= nparts(acs, :P))
println("All reporters valid: $all_valid")

# Check: all group assignments are valid
all_groups_valid = all(1 .<= subpart(acs, :part_group) .<= nparts(acs, :G)) &&
                   all(1 .<= subpart(acs, :cnt_group) .<= nparts(acs, :G))
println("All group assignments valid: $all_groups_valid")
```

    All reporters valid: true
    All group assignments valid: true

### ContactMatrixSchema

The matrix schema represents the output of the functor:

- **Objects**: G (Age Groups), E (Matrix Entries)
- **Morphisms**: `row_group: E → G`, `col_group: E → G`
- **Attributes**: `gname: G → String`, `pop: G → Float64`,
  `value: E → Float64`

``` julia
# Compute a contact matrix and convert to ACSet
pop = [6987.0, 11537.0, 35854.0, 16000.0, 9492.0]
cm = compute_matrix(survey, partition; population=pop)
cm_sym = symmetrise(cm)

acs_cm = LabelledContactMatrix(cm_sym)
println("\nContact Matrix ACSet:")
println("  Groups (G): $(nparts(acs_cm, :G))")
println("  Entries (E): $(nparts(acs_cm, :E))")
println("  Group names: $(subpart(acs_cm, :gname))")
```


    Contact Matrix ACSet:
      Groups (G): 5
      Entries (E): 25
      Group names: ["[0,5)", "[5,18)", "[18,45)", "[45,65)", "65+"]

## Functorial Data Migration

The power of ACSets is that morphisms between schemas induce **data
migration functors**. For contact data, the key morphism is coarsening:
a surjective map $f: G_{\mathrm{fine}} \to G_{\mathrm{coarse}}$ between
age-group objects.

The **left Kan extension** $\Sigma_f$ pushes data forward along $f$:

``` julia
# Start with fine partition
fine = AgePartition([0, 5, 18, 45, 65])
acs_fine = ContactSurveyACSet(survey, fine)
println("Fine partition: $(nparts(acs_fine, :G)) groups")

# Define coarsening map
coarse = AgePartition([0, 18, 65])
f = AgeMap(fine, coarse)

# Apply functorial data migration (Σ_f)
acs_coarse = migrate_coarsen(acs_fine, f)
println("Coarse partition: $(nparts(acs_coarse, :G)) groups")
println("Participants preserved: $(nparts(acs_coarse, :P))")
println("Contacts preserved: $(nparts(acs_coarse, :C))")
```

    Fine partition: 5 groups
    Coarse partition: 3 groups
    Participants preserved: 1012
    Contacts preserved: 11873

Key property: data migration preserves the number of data points
(participants and contacts) while changing only the group assignments.
No information is lost — only resolution.

``` julia
# The group assignments have been remapped
fine_groups = subpart(acs_fine, :part_group)
coarse_groups = subpart(acs_coarse, :part_group)

println("\nFirst 10 participants:")
println("  Fine groups:   $(fine_groups[1:min(10, end)])")
println("  Coarse groups: $(coarse_groups[1:min(10, end)])")
println("  (Groups $(nparts(acs_fine, :G)) → $(nparts(acs_coarse, :G)))")
```


    First 10 participants:
      Fine groups:   [2, 2, 2, 1, 1, 1, 2, 2, 1, 1]
      Coarse groups: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
      (Groups 5 → 3)

### Functoriality of Migration

Composing two coarsening maps gives the same result as the composite
map:

``` julia
# Three-level hierarchy
very_fine = AgePartition([0, 5, 10, 18, 30, 45, 65])
medium = AgePartition([0, 18, 45, 65])
very_coarse = AgePartition([0, 65])

acs_vf = ContactSurveyACSet(survey, very_fine)

# Two steps: very_fine → medium → very_coarse
f1 = AgeMap(very_fine, medium)
f2 = AgeMap(medium, very_coarse)
acs_med = migrate_coarsen(acs_vf, f1)
acs_vc_twostep = migrate_coarsen(acs_med, f2)

# One step: very_fine → very_coarse
f_composed = AgeMap(very_fine, very_coarse)
acs_vc_direct = migrate_coarsen(acs_vf, f_composed)

# Same result
println("Functoriality of data migration:")
println("  Groups match: $(nparts(acs_vc_twostep, :G) == nparts(acs_vc_direct, :G))")
println("  Part groups match: $(subpart(acs_vc_twostep, :part_group) == subpart(acs_vc_direct, :part_group))")
println("  Contact groups match: $(subpart(acs_vc_twostep, :cnt_group) == subpart(acs_vc_direct, :cnt_group))")
```

    Functoriality of data migration:
      Groups match: true
      Part groups match: true
      Contact groups match: true

## Undirected Wiring Diagrams for Composition

In the plain API, we use `⊕` for additive composition. The UWD approach
provides a **declarative** way to specify how matrices combine,
enabling:

- Named composition with explicit shared boundaries
- Automatic validation that shared boundaries are compatible
- Diagrammatic reasoning about intervention scenarios

``` julia
using Catlab.Programs: @relation

# Define a composition diagram
diagram = @relation (pop,) begin
    home(pop)
    work(pop)
    school(pop)
    other(pop)
end

println("UWD composition diagram:")
println("  Boxes: $(nparts(diagram, :Box))")
println("  Junctions: $(nparts(diagram, :Junction))")
println("  Box names: $(subpart(diagram, :name))")
```

    UWD composition diagram:
      Boxes: 4
      Junctions: 1
      Box names: [:home, :work, :school, :other]

Each box represents a contact setting. The shared junction (`:pop`)
indicates that all settings share the same population boundary.

``` julia
# Create ContactSharers for each setting
M_home = [1.8 0.3 0.4 0.1; 0.3 0.6 0.2 0.1; 0.4 0.2 1.0 0.3; 0.1 0.1 0.3 0.8]
M_work = [0.0 0.0 0.0 0.0; 0.0 0.1 1.5 0.2; 0.0 1.5 2.5 0.4; 0.0 0.2 0.4 0.3]
M_school = [3.5 2.0 0.1 0.0; 2.0 5.0 0.3 0.0; 0.1 0.3 0.2 0.0; 0.0 0.0 0.0 0.0]
M_other = [0.5 0.3 0.4 0.2; 0.3 0.5 0.8 0.3; 0.4 0.8 1.2 0.5; 0.2 0.3 0.5 0.4]

sharers = Dict(
    :home   => ContactSharer(M_home),
    :work   => ContactSharer(M_work),
    :school => ContactSharer(M_school),
    :other  => ContactSharer(M_other)
)

total = compose_uwd(diagram, sharers)
println("\nComposed matrix (via UWD):")
display(round.(total; digits=2))
```


    Composed matrix (via UWD):

    4×4 Matrix{Float64}:
     5.8  2.6  0.9  0.3
     2.6  6.2  2.8  0.6
     0.9  2.8  4.9  1.2
     0.3  0.6  1.2  1.5

### Intervention Scenarios via Diagram Modification

The UWD approach makes intervention modelling declarative — we simply
modify which boxes participate:

``` julia
# School closure: a diagram without the school box
diagram_closure = @relation (pop,) begin
    home(pop)
    work(pop)
    other(pop)
end

total_closure = compose_uwd(diagram_closure, sharers)
println("School closure scenario:")
display(round.(total_closure; digits=2))
```

    School closure scenario:

    4×4 Matrix{Float64}:
     2.3  0.6  0.8  0.3
     0.6  1.2  2.5  0.6
     0.8  2.5  4.7  1.2
     0.3  0.6  1.2  1.5

``` julia
# Partial work reduction: modify the sharer, not the diagram
sharers_lockdown = copy(sharers)
sharers_lockdown[:work] = ContactSharer(0.3 .* M_work)

total_lockdown = compose_uwd(diagram, sharers_lockdown)
println("\nLockdown scenario (70% work reduction):")
display(round.(total_lockdown; digits=2))
```


    Lockdown scenario (70% work reduction):

    4×4 Matrix{Float64}:
     5.8  2.6   0.9   0.3
     2.6  6.13  1.75  0.46
     0.9  1.75  3.15  0.92
     0.3  0.46  0.92  1.29

## Putting It All Together

The full categorical pipeline:

1.  **Survey → ACSet** (schema-enforced structure)
2.  **ACSet → ACSet** via functorial migration (coarsening)
3.  **ACSet → ContactMatrix** (compute functor)
4.  **ContactMatrix → ContactMatrix** via ⊕, ⊗, symmetrise, ↓

``` julia
# Full pipeline
partition_4 = AgePartition([0, 5, 18, 65])
pop_4 = [6987.0, 11537.0, 35854.0, 9492.0]

# Step 1: Survey to structured ACSet
acs_survey = ContactSurveyACSet(survey, partition_4)

# Step 2: Compute contact matrix
cm_full = compute_matrix(survey, partition_4; population=pop_4)

# Step 3: Symmetrise
cm_sym = symmetrise(cm_full)

# Step 4: Coarsen to children/adults (5 is a subset of [0,5,18,65])
cm_2grp = cm_sym ↓ AgePartition([0, 18])

# Step 5: Stratify for 2 regions
coupling_2 = [0.9 0.1; 0.1 0.9]
cm_regional = cm_2grp ⊗ coupling_2

println("Complete pipeline results:")
println("  Survey: $(nparts(acs_survey, :P)) participants, $(nparts(acs_survey, :C)) contacts")
println("  Full matrix: $(n_groups(cm_full)) groups")
println("  Symmetrised: $(n_groups(cm_sym)) groups")
println("  Coarsened: $(n_groups(cm_2grp)) groups")
println("  Regional: $(n_groups(cm_regional)) groups (2 regions × 2 ages)")
println("  Spectral radius: $(round(spectral_radius(cm_regional); digits=2))")
```

    Complete pipeline results:
      Survey: 1012 participants, 11873 contacts
      Full matrix: 4 groups
      Symmetrised: 4 groups
      Coarsened: 2 groups
      Regional: 4 groups (2 regions × 2 ages)
      Spectral radius: 11.18

## Summary

| Layer           | What it provides                   | Key operation        |
|-----------------|------------------------------------|----------------------|
| Plain structs   | Ergonomic API                      | `⊕`, `⊗`, `↓`        |
| ACSet schemas   | Type safety, referential integrity | `ContactSurveyACSet` |
| Data migration  | Functorial coarsening/restriction  | `migrate_coarsen`    |
| Wiring diagrams | Declarative composition            | `compose_uwd`        |

The categorical framework provides **formal guarantees** (functoriality,
idempotence, associativity) that are verified both in the Julia test
suite and in the companion Lean 4 proofs (`proofs/` directory).
