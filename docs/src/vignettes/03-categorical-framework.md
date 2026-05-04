# Categorical Framework

## Overview

This vignette dives into the **deeper categorical machinery** of ContACT.jl,
powered by [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl). While the
previous vignettes used the simplified API (plain structs + operators), the ACSet
layer provides:

1. **Schema-enforced structure**: Type-safe relationships between participants,
   contacts, and partition groups, preventing structural errors at construction.
2. **Functorial data migration** (``\Sigma_f``, ``\Delta_f``): Principled coarsening
   and restriction of partitioned data along morphisms.
3. **Undirected wiring diagrams (UWDs)**: Declarative specification of how
   setting-specific matrices compose via shared boundaries.

## Setup

```@example v03
using ContACT
using CSV
using DataFrames
using LinearAlgebra
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Programs: @relation
import ContACT: ⊕, ⊗, ↓, ↑, ▷, ↔, ρ
```

## ACSet Schemas

An **Attributed C-Set** (ACSet) is a functor from a schema category to **Set**.
ContACT.jl defines two schemas:

### ContactSurveySchema

The survey schema encodes the relational structure:

- **Objects**: P (Participants), C (Contacts), G (Partition Groups)
- **Morphisms**: `part_group: P → G`, `cnt_group: C → G`, `reporter: C → P`

This makes the foreign-key relationships explicit and type-checked.

```@example v03
# Load POLYMOD data
data_dir = joinpath(dirname(dirname(pathof(ContACT))), "data")
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

Convert to an ACSet representation:

```@example v03
partition = AgePartition([0, 5, 18, 45, 65])
acs = ContactSurveyACSet(survey, partition)

println("Survey ACSet:")
println("  Partition groups (G): $(nparts(acs, :G))")
println("  Participants (P): $(nparts(acs, :P))")
println("  Contacts (C): $(nparts(acs, :C))")
```

The ACSet enforces referential integrity: every contact's `:reporter` points to a
valid participant, and every `:part_group`/`:cnt_group` points to a valid partition group.

```@example v03
# Check: all reporters are valid participants
all_valid = all(1 .<= subpart(acs, :reporter) .<= nparts(acs, :P))
println("All reporters valid: $all_valid")

# Check: all group assignments are valid
all_groups_valid = all(1 .<= subpart(acs, :part_group) .<= nparts(acs, :G)) &&
                   all(1 .<= subpart(acs, :cnt_group) .<= nparts(acs, :G))
println("All group assignments valid: $all_groups_valid")
```

### ContactMatrixSchema

The matrix schema represents the output of the functor:

- **Objects**: G (Partition Groups), E (Matrix Entries)
- **Morphisms**: `row_group: E → G`, `col_group: E → G`
- **Attributes**: `gname: G → String`, `pop: G → Float64`, `value: E → Float64`

```@example v03
# Compute a contact matrix and convert to ACSet
pop = [6987.0, 11537.0, 35854.0, 16000.0, 9492.0]
cm = compute_matrix(survey, partition; population=pop)
cm_sym = ↔(cm)

acs_cm = LabelledContactMatrix(cm_sym)
println("\nContact Matrix ACSet:")
println("  Groups (G): $(nparts(acs_cm, :G))")
println("  Entries (E): $(nparts(acs_cm, :E))")
println("  Group names: $(subpart(acs_cm, :gname))")
```

## Functorial Data Migration

The power of ACSets is that morphisms between schemas induce **data migration
functors**. For contact data, the key morphism is coarsening: a surjective map
``f: G_{\mathrm{fine}} \to G_{\mathrm{coarse}}`` between partition-group objects.

The **left Kan extension** ``\Sigma_f`` pushes data forward along ``f``:

```@example v03
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

Key property: data migration preserves the number of data points (participants and
contacts) while changing only the group assignments. No information is lost — only
resolution.

```@example v03
# The group assignments have been remapped
fine_groups = subpart(acs_fine, :part_group)
coarse_groups = subpart(acs_coarse, :part_group)

println("\nFirst 10 participants:")
println("  Fine groups:   $(fine_groups[1:min(10, end)])")
println("  Coarse groups: $(coarse_groups[1:min(10, end)])")
println("  (Groups ``(nparts(acs_fine, :G)) → ``(nparts(acs_coarse, :G)))")
```

### Functoriality of Migration

Composing two coarsening maps gives the same result as the composite map.
This is the key property expressed with the `∘` operator:

```@example v03
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

# One step using composed map: f2 ∘ f1
f_composed = f2 ∘ f1   # type \circ<TAB>
acs_vc_direct = migrate_coarsen(acs_vf, f_composed)

# Same result — functoriality!
println("Functoriality of data migration (Σ_{g∘f} = Σ_g ∘ Σ_f):")
println("  Groups match: $(nparts(acs_vc_twostep, :G) == nparts(acs_vc_direct, :G))")
println("  Part groups match: $(subpart(acs_vc_twostep, :part_group) == subpart(acs_vc_direct, :part_group))")
println("  Contact groups match: $(subpart(acs_vc_twostep, :cnt_group) == subpart(acs_vc_direct, :cnt_group))")
```

## Undirected Wiring Diagrams for Composition

In the plain API, we use `⊕` for additive composition. The UWD approach provides
a **declarative** way to specify how matrices combine, enabling:

- Named composition with explicit shared boundaries
- Automatic validation that shared boundaries are compatible
- Diagrammatic reasoning about intervention scenarios

```@example v03
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

Each box represents a contact setting. The shared junction (`:pop`) indicates
that all settings share the same population boundary.

```@example v03
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

### Intervention Scenarios via Diagram Modification

The UWD approach makes intervention modelling declarative — we simply modify which
boxes participate:

```@example v03
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

```@example v03
# Partial work reduction: modify the sharer, not the diagram
sharers_lockdown = copy(sharers)
sharers_lockdown[:work] = ContactSharer(0.3 .* M_work)

total_lockdown = compose_uwd(diagram, sharers_lockdown)
println("\nLockdown scenario (70% work reduction):")
display(round.(total_lockdown; digits=2))
```

## Putting It All Together

The full categorical pipeline expressed with Unicode operators:

1. **Survey → ACSet** (schema-enforced structure)
2. **ACSet → ACSet** via functorial migration (coarsening)
3. **Survey ▷ Partition** (compute functor)
4. **ContactMatrix ↓ ∘ ⊕ ⊗** (algebraic manipulation)

```@example v03
# Full pipeline using operators
partition_4 = AgePartition([0, 5, 18, 65])
pop_4 = [6987.0, 11537.0, 35854.0, 9492.0]

# Step 1: Survey to structured ACSet
acs_survey = ContactSurveyACSet(survey, partition_4)

# Step 2: Apply the functor (▷)
cm_full = compute_matrix(survey, partition_4; population=pop_4)

# Step 3: Symmetrise via ↔ (idempotent reciprocity projection)
cm_sym = ↔(cm_full)

# Step 4: Coarsen via ↓
cm_2grp = cm_sym ↓ AgePartition([0, 18])

# Step 5: Stratify via ⊗
coupling_2 = [0.9 0.1; 0.1 0.9]
cm_regional = cm_2grp ⊗ coupling_2

println("Complete pipeline results:")
println("  Survey: ``(nparts(acs_survey, :P)) participants, ``(nparts(acs_survey, :C)) contacts")
println("  Full matrix: ``(n_groups(cm_full)) groups, ρ = ``(round(ρ(cm_full); digits=2))")
println("  Symmetrised: ``(n_groups(cm_sym)) groups, ρ = ``(round(ρ(cm_sym); digits=2))")
println("  Coarsened: ``(n_groups(cm_2grp)) groups, ρ = ``(round(ρ(cm_2grp); digits=2))")
println("  Regional: ``(n_groups(cm_regional)) groups, ρ = ``(round(ρ(cm_regional); digits=2))")
```

## Summary

| Operator | LaTeX | Category-theoretic role |
|----------|-------|------------------------|
| `⊕` | `\\oplus` | Monoidal product (commutative monoid) |
| `⊗` | `\\otimes` | Kronecker/stratification functor |
| `↓` | `\\downarrow` | Left Kan extension (coarsening) |
| `↑` | `\\uparrow` | Parameterised refinement (with prior) |
| `▷` | `\\triangleright` | Functor application (survey → matrix) |
| `∘` | `\\circ` | Morphism composition (PartitionMap; `AgeMap` is the age alias) |
| `↔` | `\\leftrightarrow` | Reciprocity projection (symmetrisation) |
| `ρ` | `\\rho` | Spectral radius (R₀ proxy) |

The categorical framework provides **formal guarantees** (functoriality, idempotence,
associativity) that are verified both in the Julia test suite and in the companion
Lean 4 proofs (`proofs/` directory).
