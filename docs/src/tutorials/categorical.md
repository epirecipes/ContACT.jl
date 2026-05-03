# Categorical Framework

This tutorial demonstrates the deeper categorical machinery powered by
[Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl): ACSet schemas,
functorial data migration, and undirected wiring diagrams.

## ACSet Schemas

An **Attributed C-Set** (ACSet) is a functor from a schema category to **Set**.
ContACT.jl defines schemas that enforce structural invariants at construction time.

### Survey Schema

The survey schema captures participants, contacts, and partition groups with
explicit morphisms:

```@example categorical
using ContACT
using DataFrames
using CSV
using Catlab
using Catlab.CategoricalAlgebra
import ContACT: ⊕, ⊗, ↓

# Load POLYMOD data
data_dir = joinpath(dirname(dirname(pathof(ContACT))), "data")
participants = CSV.read(joinpath(data_dir, "polymod_uk_participants.csv"), DataFrame)
contacts = CSV.read(joinpath(data_dir, "polymod_uk_contacts.csv"), DataFrame)

parse_age(x::AbstractString) = x == "NA" ? missing : parse(Float64, x)
parse_age(x::Real) = Float64(x)
parse_age(::Missing) = missing

rename!(participants, :part_age_exact => :part_age)
participants.part_age = Float64.(participants.part_age)
cnt_exact = parse_age.(contacts.cnt_age_exact)
cnt_min = parse_age.(contacts.cnt_age_est_min)
cnt_max = parse_age.(contacts.cnt_age_est_max)
contacts.cnt_age = coalesce.(cnt_exact, (cnt_min .+ cnt_max) ./ 2)
select!(contacts, [:part_id, :cnt_age])
dropmissing!(contacts, :cnt_age)

survey = ContactSurvey(participants, contacts)

# Convert to ACSet
partition = AgePartition([0, 5, 18, 45, 65])
acs = ContactSurveyACSet(survey, partition)

println("Survey ACSet:")
println("  Partition groups (G): $(nparts(acs, :G))")
println("  Participants (P): $(nparts(acs, :P))")
println("  Contacts (C): $(nparts(acs, :C))")
nothing # hide
```

The ACSet enforces referential integrity:

```@example categorical
all_reporters_valid = all(1 .<= subpart(acs, :reporter) .<= nparts(acs, :P))
all_groups_valid = all(1 .<= subpart(acs, :part_group) .<= nparts(acs, :G)) &&
                   all(1 .<= subpart(acs, :cnt_group) .<= nparts(acs, :G))
println("All reporters valid: $all_reporters_valid")
println("All group assignments valid: $all_groups_valid")
nothing # hide
```

### Contact Matrix Schema

```@example categorical
pop = [6987.0, 11537.0, 35854.0, 16000.0, 9492.0]
cm = compute_matrix(survey, partition; population=pop)
cm_sym = ↔(cm)

acs_cm = LabelledContactMatrix(cm_sym)
println("Contact Matrix ACSet:")
println("  Groups (G): $(nparts(acs_cm, :G))")
println("  Entries (E): $(nparts(acs_cm, :E))")
println("  Labels: $(subpart(acs_cm, :gname))")
nothing # hide
```

## Functorial Data Migration

A surjective map ``f: G_{\text{fine}} \to G_{\text{coarse}}`` between partition-group
objects induces a **left Kan extension** ``\Sigma_f`` that pushes survey data
forward:

```@example categorical
fine = AgePartition([0, 5, 18, 45, 65])
acs_fine = ContactSurveyACSet(survey, fine)

coarse = AgePartition([0, 18, 65])
f = AgeMap(fine, coarse)
acs_coarse = migrate_coarsen(acs_fine, f)

println("Migration Σ_f:")
println("  Groups: $(nparts(acs_fine, :G)) → $(nparts(acs_coarse, :G))")
println("  Participants preserved: $(nparts(acs_coarse, :P))")
println("  Contacts preserved: $(nparts(acs_coarse, :C))")
nothing # hide
```

### Functoriality

Composing two migrations equals a single migration along the composite map:

```@example categorical
very_fine = AgePartition([0, 5, 10, 18, 30, 45, 65])
medium = AgePartition([0, 18, 45, 65])
very_coarse = AgePartition([0, 65])

acs_vf = ContactSurveyACSet(survey, very_fine)

# Two-step migration
f1 = AgeMap(very_fine, medium)
f2 = AgeMap(medium, very_coarse)
acs_twostep = migrate_coarsen(migrate_coarsen(acs_vf, f1), f2)

# One-step migration
f_composed = AgeMap(very_fine, very_coarse)
acs_direct = migrate_coarsen(acs_vf, f_composed)

println("Functoriality: Σ_{f₂} ∘ Σ_{f₁} = Σ_{f₂∘f₁}")
println("  Part groups match: $(subpart(acs_twostep, :part_group) == subpart(acs_direct, :part_group))")
println("  Contact groups match: $(subpart(acs_twostep, :cnt_group) == subpart(acs_direct, :cnt_group))")
nothing # hide
```

## Undirected Wiring Diagrams

UWDs provide a **declarative** specification of how setting-specific matrices
compose. Each box is a setting; shared junctions represent compatible boundaries.

```@example categorical
using Catlab.Programs: @relation

diagram = @relation (pop,) begin
    home(pop)
    work(pop)
    school(pop)
    other(pop)
end

println("UWD: $(nparts(diagram, :Box)) boxes, $(nparts(diagram, :Junction)) junctions")
println("Box names: $(subpart(diagram, :name))")
nothing # hide
```

```@example categorical
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
round.(total; digits=2)
```

### Intervention via Diagram Modification

```@example categorical
# School closure: simply omit the box
diagram_closure = @relation (pop,) begin
    home(pop)
    work(pop)
    other(pop)
end

total_closure = compose_uwd(diagram_closure, sharers)
println("School closure matrix:")
round.(total_closure; digits=2)
```

## Full Pipeline

```@example categorical
partition_4 = AgePartition([0, 5, 18, 65])
pop_4 = [6987.0, 11537.0, 35854.0, 9492.0]

acs_survey = ContactSurveyACSet(survey, partition_4)
cm_full = compute_matrix(survey, partition_4; population=pop_4)
cm_sym2 = ↔(cm_full)
cm_2grp = cm_sym2 ↓ AgePartition([0, 18])
cm_regional = cm_2grp ⊗ [0.9 0.1; 0.1 0.9]

println("Pipeline: Survey → ACSet → Matrix → ↔ → ↓ → ⊗")
println("  Input: $(nparts(acs_survey, :P)) participants")
println("  Final: $(n_groups(cm_regional)) groups")
println("  Spectral radius: $(round(ρ(cm_regional); digits=2))")
nothing # hide
```
