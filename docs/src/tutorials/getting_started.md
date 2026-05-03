# Getting Started

This tutorial introduces ContACT.jl using the POLYMOD study — a multi-country
social contact survey that is the gold standard for parameterising infectious
disease models.

## Loading Survey Data

ContACT.jl works with `ContactSurvey` objects containing participant and contact
DataFrames. The survey itself only requires a participant identifier; the
partition you choose specifies which participant/contact columns define groups.

```@example tutorial
using ContACT
using CSV
using DataFrames

# Load POLYMOD UK data bundled with the package
data_dir = joinpath(dirname(dirname(pathof(ContACT))), "data")
participants = CSV.read(joinpath(data_dir, "polymod_uk_participants.csv"), DataFrame)
contacts = CSV.read(joinpath(data_dir, "polymod_uk_contacts.csv"), DataFrame)

println("Participants: $(nrow(participants))")
println("Contacts: $(nrow(contacts))")
nothing # hide
```

## Preparing the Survey

POLYMOD stores contact ages as strings (with "NA" for missing). We parse and
harmonise:

```@example tutorial
# Parse string ages to Float64
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
```

## Computing a Contact Matrix

The central operation is the **restricted functor** from surveys to contact
matrices. We fix a partition to make the operation deterministic; age is the
default interval-valued partition:

```@example tutorial
# Standard 5-year age bands
partition = AgePartition(
    [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
)

# UK 2005 population (thousands, approximate)
uk_pop = [3433.0, 3554.0, 3824.0, 3956.0, 3757.0, 3520.0, 4009.0,
           4405.0, 4548.0, 4187.0, 3883.0, 3589.0, 3117.0, 2765.0,
           2425.0, 3752.0]

# Apply the survey→matrix functor using ▷ (type \triangleright<TAB>)
cm = survey ▷ partition
println("Contact matrix: $(n_groups(cm)) age groups")
println("ρ(M) ∝ R₀: $(round(ρ(cm); digits=2))")
nothing # hide
```

The result bundles the matrix with its partition, population, and semantics:

```@example tutorial
println("Semantics: $(cm.semantics)")
println("Age groups: $(age_labels(cm))")
nothing # hide
```

## Non-age partitions

The same functor works for categorical survey variables such as sex, region,
occupation, or products of those variables. The partition carries the semantic
dimension (`:sex`) and the column names used in the participant and contact
tables:

```@example tutorial
participants_sex = DataFrame(part_id=[1, 2, 3], part_sex=["F", "M", "F"])
contacts_sex = DataFrame(
    part_id=[1, 1, 2, 3],
    cnt_sex=["M", "F", "F", "F"],
)
survey_sex = ContactSurvey(participants_sex, contacts_sex)

sex = CategoricalPartition(:sex;
    participant_col=:part_sex,
    contact_col=:cnt_sex,
    levels=["F", "M"],
    labels=["female", "male"],
)

cm_sex = survey_sex ▷ sex
println(group_labels(cm_sex))
matrix(cm_sex)
```

Product partitions combine grouping variables without inventing artificial
numeric labels:

```@example tutorial
region = CategoricalPartition(:region;
    participant_col=:part_region,
    contact_col=:cnt_region,
    levels=["north", "south"],
)
sex_region = sex × region
group_labels(sex_region)
```

## Symmetrisation

Raw survey matrices are asymmetric. Symmetrisation enforces reciprocity
(``M_{ij} N_j = M_{ji} N_i``) and is **idempotent**. Use `↔`
(`\leftrightarrow<TAB>`) as the reciprocity projection:

```@example tutorial
cm_sym = ↔(cm)

# Check idempotence
cm_sym2 = ↔(cm_sym)
println("Idempotent: $(matrix(cm_sym2) ≈ matrix(cm_sym))")
nothing # hide
```

```@example tutorial
# Verify reciprocity
M = matrix(cm_sym)
N = population(cm_sym)
max_violation = maximum(
    abs(M[i,j] * N[j] - M[j,i] * N[i])
    for i in 1:n_groups(cm_sym), j in 1:n_groups(cm_sym)
)
println("Max reciprocity violation: $(max_violation)")
nothing # hide
```

## Coarsening

Coarsening is the left Kan extension along a surjective partition map. For
interval partitions such as age, the coarse limits must be a subset of the fine
limits:

```@example tutorial
coarse = AgePartition([0, 15, 45, 65])
cm_coarse = cm ↓ coarse

println("Coarsened: $(n_groups(cm_coarse)) groups")
println("Labels: $(age_labels(cm_coarse))")
round.(matrix(cm_coarse); digits=2)
```

### Functoriality

Coarsening in two steps gives the same result as one step — the defining
property verified in our Lean proofs. Express with `∘` (map composition):

```@example tutorial
mid = AgePartition([0, 15, 45, 65])
final_part = AgePartition([0, 45])

# Compose maps with ∘ (type \circ<TAB>)
f = AgeMap(partition, mid)
g = AgeMap(mid, final_part)
h = g ∘ f

via_mid = (cm ↓ f) ↓ g
direct = cm ↓ h

println("cm ↓ (g ∘ f) ≈ (cm ↓ f) ↓ g: $(matrix(via_mid) ≈ matrix(direct))")
round.(matrix(direct); digits=2)
```
