# Categorical Analysis of the Reconnect UK Contact Survey

## Overview

This vignette demonstrates a **categorical analysis** of the
[Reconnect UK social contact survey](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1005038)
(Goodfellow et al., 2026) using ContACT.jl. The Reconnect study collected
individual-level contact data from 13,238 participants across the UK,
recording details of 125,671 social contacts including age, setting (home,
work, school, other), and socio-demographic characteristics.

We show how each step of the analysis corresponds to a morphism or functor in
the categorical framework:

| Step                        | Categorical Operation                          | Operator |
|:----------------------------|:-----------------------------------------------|:---------|
| Survey → contact matrix     | Functor `F: Survey(P) → Contact`               | `▷`      |
| Symmetrisation              | Idempotent endomorphism                        | `↔`      |
| Setting composition         | Commutative monoid `⊕`                         | `⊕`      |
| Age coarsening              | Left Kan extension along `FinFunction`          | `↓`      |
| NGM construction            | Diagonal scaling by population                 | —        |

!!! note "Methodological Note"
    This vignette computes **simple day-weighted arithmetic mean contacts**
    per participant group. The published Reconnect analysis uses a more
    sophisticated negative binomial MLE bootstrap estimator with
    post-stratification weights and large-group contact allocation.
    Results here reflect only individually reported contacts with known
    contact age (~50,665 of 125,671 contacts).

## Setup

```@example v08
using ContACT
using CSV
using DataFrames
using LinearAlgebra
using Downloads: download
import ContACT: ×
```

## 1. Loading the Reconnect Data

The Reconnect survey data are available from
[Zenodo (DOI: 10.5281/zenodo.17339866)](https://zenodo.org/records/17339866).

```@example v08
zenodo_base = "https://zenodo.org/records/17339866/files"
data_dir = mktempdir()

filenames = [
    "reconnect_participant_common.csv",
    "reconnect_participant_extra.csv",
    "reconnect_contact_common.csv",
    "reconnect_contact_extra.csv",
    "reconnect_sday.csv",
]

for f in filenames
    dest = joinpath(data_dir, f)
    isfile(dest) || download("$zenodo_base/$f", dest)
end
println("Data downloaded to: $data_dir")
```

```@example v08
part_common = CSV.read(
    joinpath(data_dir, "reconnect_participant_common.csv"), DataFrame;
    missingstring = "NA",
)
part_extra = CSV.read(
    joinpath(data_dir, "reconnect_participant_extra.csv"), DataFrame;
    missingstring = "NA",
)
sday = CSV.read(joinpath(data_dir, "reconnect_sday.csv"), DataFrame;
    missingstring = "NA",
)
cnt_common = CSV.read(
    joinpath(data_dir, "reconnect_contact_common.csv"), DataFrame;
    missingstring = "NA",
)
cnt_extra = CSV.read(
    joinpath(data_dir, "reconnect_contact_extra.csv"), DataFrame;
    missingstring = "NA",
)

println("Participants (common): $(nrow(part_common))")
println("Contacts (common):     $(nrow(cnt_common))")
```

## 2. Data Preparation

Join participant tables and compute day-of-week weights.

```@example v08
participants = innerjoin(part_common, part_extra, on = :part_id)
participants = innerjoin(participants, sday, on = :part_id)

# Drop any participant with missing day-of-week
dropmissing!(participants, :dayofweek)

# Day type: dayofweek 0=Sunday, 1=Monday, ..., 6=Saturday
# Weekend = 0 (Sunday) and 6 (Saturday)
participants.day_type = ifelse.(
    participants.dayofweek .∈ Ref([0, 6]),
    "weekend", "weekday"
)

# Day weights: reweight to 5/7 weekday, 2/7 weekend
day_counts = combine(groupby(participants, :day_type), nrow => :n)
total_n = nrow(participants)
day_counts.sample_prop = day_counts.n ./ total_n
day_counts.target_prop = ifelse.(
    day_counts.day_type .== "weekday", 5 / 7, 2 / 7
)
day_counts.day_weight = day_counts.target_prop ./ day_counts.sample_prop

participants = innerjoin(
    participants,
    select(day_counts, :day_type, :day_weight),
    on = :day_type,
)
println("Day weight distribution:")
println(select(day_counts, :day_type, :n, :sample_prop, :day_weight))
```

Prepare the contact data.

```@example v08
contacts = innerjoin(cnt_common, cnt_extra, on = [:cont_id, :part_id])

# Rename age columns to match ContACT conventions
rename!(participants, :part_age_exact => :part_age)
participants.part_age = Float64.(participants.part_age)
rename!(contacts, :cnt_age_exact => :cnt_age)

# Retain only contacts with known exact age
n_before = nrow(contacts)
dropmissing!(contacts, :cnt_age)
contacts.cnt_age = Float64.(collect(contacts.cnt_age))

println(
    "Contacts with known age: $(nrow(contacts)) of $n_before " *
    "($(round(100 * nrow(contacts) / n_before, digits=1))%)"
)
```

```@example v08
# Keep only contacts whose participants are in our dataset
valid_ids = Set(participants.part_id)
filter!(row -> row.part_id in valid_ids, contacts)

# Build the ContactSurvey
survey = ContactSurvey(participants, contacts)
println(
    "Survey: $(nrow(survey.participants)) participants, " *
    "$(nrow(survey.contacts)) contacts"
)
```

## 3. Age-Stratified Contact Matrix

We apply the **survey-to-matrix functor** with a 16-group age partition
matching the Reconnect study's published age bands.

```@example v08
age_breaks = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
age_part = AgePartition(age_breaks)

# The functor maps survey data to a ContactMatrix, weighted by day_weight
cm_age = compute_matrix(survey, age_part; weights = :day_weight)

println("Contact matrix: $(n_groups(cm_age)) × $(n_groups(cm_age))")
println(
    "\nMean contacts with known age, by participant age group (column sums):"
)
M = matrix(cm_age)
for (j, lab) in enumerate(group_labels(cm_age))
    println("  $lab: $(round(sum(M[:, j]), digits=2))")
end
```

!!! note
    These column sums represent contacts with **known contact age** only
    (~50,665 of 125,671 total contacts). The paper reports higher mean
    daily contacts (9.1) because it includes large-group contacts and
    contacts with imputed ages.

```@example v08
# Display as a compact rounded matrix
println("Contact matrix (rows=contact age, cols=participant age):")
show(
    IOContext(stdout, :compact => true, :limit => true),
    "text/plain",
    round.(matrix(cm_age), digits = 3),
)
println()
```

## 4. Symmetrisation

The **symmetrisation morphism** `↔` enforces the reciprocity constraint.

```math
M^{\text{sym}}_{ij} = \frac{M_{ij} \cdot N_j + M_{ji} \cdot N_i}{2\,N_j}
```

This is an **idempotent** endomorphism: `↔(↔(M)) = ↔(M)`.

```@example v08
cm_sym = ↔(cm_age)

# Verify reciprocity: M[i,j]*N[j] ≈ M[j,i]*N[i]
M_s = matrix(cm_sym)
N = population(cm_sym)
max_reciprocity_error = maximum(
    abs(M_s[i, j] * N[j] - M_s[j, i] * N[i]) for i in 1:16 for j in 1:16
)
println(
    "Maximum reciprocity violation: " *
    "$(round(max_reciprocity_error, sigdigits=3))"
)

# Verify idempotence
cm_sym2 = ↔(cm_sym)
println(
    "Idempotence check: max|↔(↔(M)) − ↔(M)| = " *
    "$(maximum(abs.(matrix(cm_sym2) .- M_s)))"
)
```

## 5. Setting-Specific Composition

The total contact matrix is the **additive composition** of setting-specific
matrices: `M_total = M_home ⊕ M_work ⊕ M_school ⊕ M_other`.

```@example v08
# Verify that location categories partition the contact set
location_counts = combine(
    groupby(DataFrame(cnt_location = contacts.cnt_location), :cnt_location),
    nrow => :n,
)
println("Contact locations:")
println(location_counts)
```

```@example v08
settings = ["Home", "Work", "School", "Other"]
cms = Dict{String,ContactMatrix}()

for setting in settings
    cnt_setting = filter(row -> row.cnt_location == setting, contacts)
    survey_setting = ContactSurvey(participants, cnt_setting)
    cms[setting] = compute_matrix(
        survey_setting, age_part; weights = :day_weight
    )
end

# Additive composition: ⊕
cm_composed = cms["Home"] ⊕ cms["Work"] ⊕ cms["School"] ⊕ cms["Other"]

# Verify the composition property
composition_error = maximum(abs.(matrix(cm_composed) .- matrix(cm_age)))
println(
    "Composition check: max|Home ⊕ Work ⊕ School ⊕ Other − Total| = " *
    "$(round(composition_error, sigdigits=3))"
)
```

```@example v08
println("Mean known-age contacts per person by setting:")
for setting in settings
    total = round(sum(matrix(cms[setting])), digits = 2)
    println("  $setting: $total")
end
total_all = round(sum(matrix(cm_age)), digits = 2)
println("  Total:  $total_all")
```

## 6. Ethnicity-Stratified Matrix

We use a **categorical partition** with complete-case analysis: only
participants and contacts with ethnicity in the five main categories
are included.

```@example v08
eth_levels = ["Asian", "Black", "Mixed", "Other", "White"]

eth_part = CategoricalPartition(
    :ethnicity,
    eth_levels;
    participant_col = :part_ethnicity,
    contact_col = :cnt_ethnicity,
)

# Complete case: filter participants and contacts
participants_eth = filter(
    row ->
        !ismissing(row.part_ethnicity) && row.part_ethnicity in eth_levels,
    participants,
)

valid_eth_ids = Set(participants_eth.part_id)
contacts_eth = filter(
    row ->
        row.part_id in valid_eth_ids &&
            !ismissing(row.cnt_ethnicity) &&
            row.cnt_ethnicity in eth_levels,
    contacts,
)

survey_eth = ContactSurvey(participants_eth, contacts_eth)
cm_eth = compute_matrix(survey_eth, eth_part; weights = :day_weight)

println(
    "Ethnicity matrix ($(nrow(participants_eth)) participants, " *
    "$(nrow(contacts_eth)) contacts):"
)
```

```@example v08
M_eth = matrix(cm_eth)
labels_eth = group_labels(cm_eth)
print(lpad("", 8))
for lab in labels_eth
    print(rpad(lab, 10))
end
println()
for (j, lab) in enumerate(labels_eth)
    print(rpad(lab, 8))
    for i in eachindex(labels_eth)
        print(rpad(round(M_eth[i, j], digits = 3), 10))
    end
    println()
end
```

## 7. Coarsening: Fine → Coarse Age Groups

Coarsening is a left Kan extension along a `FinFunction` that maps fine
groups to coarse groups.

!!! warning "Partition Compatibility"
    The 16-group partition `[0, 5, 10, 15, 20, ...]` does **not** refine
    `[0, 18, 65]` because age 18 falls inside the `[15, 20)` bin.

### Direct computation (exact)

```@example v08
coarse_part = AgePartition([0, 18, 65])
cm_coarse_direct = compute_matrix(survey, coarse_part; weights = :day_weight)

println("Direct 3×3 age matrix (0-17, 18-64, 65+):")
println(round.(matrix(cm_coarse_direct), digits = 4))
```

### Approximate coarsening (from 16 groups)

```@example v08
# Bin [15,20) straddles the 18 boundary → assigned to children (majority)
assignments = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3]
age_map = PartitionMap(age_part, coarse_part, assignments)

cm_coarse_approx = cm_age ↓ age_map

println("Approximate 3×3 (15-19 → children):")
println(round.(matrix(cm_coarse_approx), digits = 4))

println("\nDifference (approx − direct):")
println(round.(matrix(cm_coarse_approx) .- matrix(cm_coarse_direct), digits = 4))
```

## 8. Next-Generation Matrix and R₀

```@example v08
τ = 0.05
K = next_generation_matrix(cm_sym; transmissibility = τ)
R0_val = R₀(cm_sym; transmissibility = τ)

println("R₀ = $(round(R0_val, digits=3)) (with τ = $τ)")
println("Spectral radius via ρ(cm): $(round(ρ(cm_sym) * τ, digits=3))")
```

```@example v08
pop_fracs = population(cm_sym) ./ sum(population(cm_sym))
τ_vec = solve_final_size_vector(K, pop_fracs)

println("Per-group final attack fractions (τ = $τ):")
for (i, lab) in enumerate(group_labels(cm_sym))
    println("  $lab: $(round(τ_vec[i] * 100, digits=1))%")
end
println(
    "Overall attack fraction: " *
    "$(round(sum(τ_vec .* pop_fracs) * 100, digits=1))%"
)
```

## 9. Epidemic Bounds from Partial Information

```@example v08
# Bounds on R₀ from only the row/column marginals of the NGM
bounds = r0_bounds(K)
println("R₀ bounds from marginals alone:")
println("  Lower bound: $(round(bounds.lower, digits=3))")
println("  Upper bound: $(round(bounds.upper, digits=3))")
println("  True R₀:     $(round(R0_val, digits=3))")
```

## 10. Summary

| Section | Operation | Key Result |
|:--------|:----------|:-----------|
| §3 | Functor `▷` | Survey → 16×16 age contact matrix |
| §4 | Symmetrisation `↔` | Reciprocity enforced; idempotent |
| §5 | Composition `⊕` | Home ⊕ Work ⊕ School ⊕ Other = Total (exact) |
| §6 | Categorical partition | 5×5 ethnicity matrix; assortativity visible |
| §7 | Coarsening `↓` | Requires compatible partitions |
| §8 | NGM + R₀ | Spectral radius of τC |
| §9 | Epidemic bounds | R₀ bracketed from marginals alone |

## References

- Goodfellow I, et al. (2026). "Social contact patterns in the UK from the
  Reconnect study." *PLOS Medicine*.
- Mossong J, et al. (2008). "Social contacts and mixing patterns relevant to
  the spread of infectious diseases." *PLOS Medicine*, 5(3), e74.
- Britton T, et al. (2025). "Bounds on R₀ and final epidemic size."
  arXiv:2602.23885v2.
