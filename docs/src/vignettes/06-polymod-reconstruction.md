# Real-Data Reconstruction with POLYMOD

## Overview

This vignette demonstrates the **constrained generalized lift** pipeline using
real contact survey data from the POLYMOD UK study (Mossong et al., 2008). We
reconstruct a fully stratified contact matrix over `age × daytype` from partial
observations where contact ages are known but the day-type of the contact event
is not directly available for the contact side.

This illustrates the core CoMix/SEP workflow:
1. Observe a rectangular source-stratified matrix (contacts classified by target
   age only, but participants classified by age × day-type)
2. Use a reciprocal base age-only matrix as constraint
3. Reconstruct the full square matrix over the product partition
4. Explore the family of valid reconstructions via q-parameters and MCMC

## Setup

```@example v06
using ContACT
using CSV
using DataFrames
using LinearAlgebra
using Random
import ContACT: ×
```

## Loading POLYMOD UK Data

```@example v06
data_dir = joinpath(dirname(dirname(pathof(ContACT))), "data")
participants = CSV.read(joinpath(data_dir, "polymod_uk_participants.csv"), DataFrame)
contacts = CSV.read(joinpath(data_dir, "polymod_uk_contacts.csv"), DataFrame)

println("Participants: ``(nrow(participants)), Contacts: ``(nrow(contacts))")
```

## Data Preparation

We create a `daytype` covariate: "weekday" (Mon-Fri) vs "weekend" (Sat-Sun).
The POLYMOD `dayofweek` column uses 0=Sunday through 6=Saturday.

```@example v06
# Parse ages
parse_age(x::AbstractString) = x == "NA" ? missing : parse(Float64, x)
parse_age(x::Real) = Float64(x)
parse_age(::Missing) = missing

# Participant ages
rename!(participants, :part_age_exact => :part_age)
participants.part_age = Float64.(participants.part_age)

# Day type: 0=Sun, 6=Sat -> weekend; 1-5 -> weekday
function classify_daytype(d)
    d in ("0", "6") && return "weekend"
    d == "NA" && return missing
    return "weekday"
end
participants.part_daytype = classify_daytype.(participants.dayofweek)

# Contact ages (use exact or midpoint of range)
cnt_exact = parse_age.(contacts.cnt_age_exact)
cnt_min = parse_age.(contacts.cnt_age_est_min)
cnt_max = parse_age.(contacts.cnt_age_est_max)
contacts.cnt_age = coalesce.(cnt_exact, (cnt_min .+ cnt_max) ./ 2)
select!(contacts, [:part_id, :cnt_age])
dropmissing!(contacts, :cnt_age)

# Drop participants with missing daytype
dropmissing!(participants, :part_daytype)

# Rebuild contacts to only include known participants
valid_ids = Set(participants.part_id)
filter!(row -> row.part_id in valid_ids, contacts)

survey = ContactSurvey(participants, contacts)
println("Survey: ``(nrow(survey.participants)) participants, ``(nrow(survey.contacts)) contacts")
```

## Defining Partitions

We define three age groups and a weekday/weekend partition, then form their
product. The key categorical insight: contacts are only classifiable by the
**target age** partition, but participants carry the full product partition
information.

```@example v06
# Coarse age partition: children, adults, elderly
age_partition = AgePartition([0, 18, 65])
println("Age groups: ", group_labels(age_partition))

# Day-type partition (participant-only covariate)
daytype_partition = CategoricalPartition(
    :daytype,
    ["weekday", "weekend"];
    participant_col=:part_daytype,
    contact_col=:cnt_daytype  # not available for contacts
)
println("Day-type groups: ", group_labels(daytype_partition))

# Product partition: age x daytype
product = age_partition × daytype_partition
println("Product groups ($(n_groups(product))): ", group_labels(product))
```

## Step 1: Compute the Source-Stratified Matrix

The source-stratified matrix has **target rows** indexed by the base age
partition (3 groups) and **source columns** indexed by the full product
partition (6 groups = 3 ages x 2 day-types). This is a rectangular matrix.

```@example v06
# UK 2005 approximate population in product groups
# (proportional split: ~71.4% weekday, ~28.6% weekend for each age group)
uk_pop_age = [11000.0, 33000.0, 9500.0]  # thousands: <18, 18-65, 65+
weekday_frac = 5/7
product_pop = Float64[]
for p in uk_pop_age
    push!(product_pop, p * weekday_frac)
    push!(product_pop, p * (1 - weekday_frac))
end
println("Product population: ", round.(product_pop; digits=0))

# Compute source-stratified matrix
intermediate = compute_source_stratified_matrix(
    survey, age_partition, product;
    population=product_pop
)
println("\nSource-stratified matrix ($(size(intermediate))):")
println("  Target groups (rows): ", target_group_labels(intermediate))
println("  Source groups (cols): ", source_group_labels(intermediate))
println("\nMean contacts per source group member (rows=target age):")
display(round.(matrix(intermediate); digits=2))
```

## Step 2: Compute the Base Age-Only Matrix

The base matrix provides the reciprocal constraint -- it captures the age
structure of mixing balanced by population.

```@example v06
# Compute the standard contact matrix over age only, with known population
cm_age = compute_matrix(survey, age_partition; population=uk_pop_age)

# Symmetrise to enforce reciprocity (prefix operator)
cm_age_sym = ↔(cm_age)
println("Base matrix (``(n_groups(cm_age_sym)) x ``(n_groups(cm_age_sym))):")
display(round.(matrix(cm_age_sym); digits=2))
println("\nSpectral radius rho: ", round(ρ(cm_age_sym); digits=2))
```

## Step 3: Proportionate Reconstruction (q = 0)

The constrained generalized lift reconstructs the full 6x6 matrix assuming
proportionate mixing within each age block (the default, q = 0).

```@example v06
# Build the lift specification
source_map = PartitionMap(product, age_partition)  # projection: age x daytype -> age
spec = ConstrainedGeneralizedLift(product, intermediate; source_map=source_map)

# Apply the lift (boxtimes operator: base boxtimes spec)
cm_full = cm_age_sym ⊠ spec
println("Reconstructed matrix (``(n_groups(cm_full)) x ``(n_groups(cm_full))):")
display(round.(matrix(cm_full); digits=2))
println("\nGroup labels: ", group_labels(cm_full.partition))
println("Spectral radius: ", round(ρ(cm_full); digits=2))
```

## Step 4: Coarsening Consistency Check

A key categorical property: coarsening the reconstructed full matrix back to
the age partition should recover the base matrix.

```@example v06
# Coarsen back to age (downarrow operator)
cm_coarsened = cm_full ↓ age_partition
println("Coarsened back to age-only:")
display(round.(matrix(cm_coarsened); digits=2))

println("\nOriginal base matrix:")
display(round.(matrix(cm_age_sym); digits=2))

# Check: spectral radii should match
println("\nrho(full) = ", round(ρ(cm_full); digits=4))
println("rho(base) = ", round(ρ(cm_age_sym); digits=4))
```

## Step 5: Exploring Assortativity with q-Parameters

The proportionate reconstruction (q = 0) is just one member of a family.
We can explore assortative mixing where participants preferentially contact
others with the same day-type.

```@example v06
# Global q > 0: assortative (prefer same daytype)
params_assort = BlockAssortativityParams(q=Dict(:daytype => 0.5))
pspec = ParameterizedConstrainedLift(spec; default_params=params_assort)
cm_assort = cm_age_sym ⊠ pspec

println("Assortative (q=0.5) reconstruction:")
println("  Assortativity index: ", round(assortativity_index(cm_assort, :daytype); digits=4))

# Compare with proportionate
println("  Proportionate AI:    ", round(assortativity_index(cm_full, :daytype); digits=4))

# Disassortative
params_disassort = BlockAssortativityParams(q=Dict(:daytype => -0.3))
pspec_dis = ParameterizedConstrainedLift(spec; default_params=params_disassort)
cm_disassort = cm_age_sym ⊠ pspec_dis
println("  Disassortative AI:   ", round(assortativity_index(cm_disassort, :daytype); digits=4))
```

## Step 6: Sampling the Feasible Space

Rejection sampling explores the full family of valid reconstructions over
q-parameter space.

```@example v06
Random.seed!(42)

# sample_constrained_lifts(base, spec, n; bounds) returns Vector{(params, matrix)}
samples = sample_constrained_lifts(cm_age_sym, spec, 100; bounds=(-0.8, 0.8))
println("Generated $(length(samples)) feasible reconstructions")

# Summary statistics
ais = [assortativity_index(cm, :daytype) for (_, cm) in samples]
rhos = [ρ(cm) for (_, cm) in samples]
println("Assortativity index range: [``(round(minimum(ais); digits=3)), ``(round(maximum(ais); digits=3))]")
println("Spectral radius range: [``(round(minimum(rhos); digits=3)), ``(round(maximum(rhos); digits=3))]")
```

## Step 7: MCMC with Targeted Density

For more efficient exploration, we use Metropolis-Hastings MCMC with a
log-density that favours reconstructions with moderate assortativity (e.g.
matching an observed mixing pattern).

```@example v06
# Target: moderate assortativity (AI ~ 0.3)
target_ai = 0.3

# log_density receives (ContactMatrix, block_params) -> Float64
log_density = (cm, _params) -> -50.0 * (assortativity_index(cm, :daytype) - target_ai)^2

# Run MCMC (auto-detects dimensions and builds QParameterSpace internally)
result = mcmc_constrained_lifts(
    cm_age_sym, spec, 500;
    bounds=(-0.8, 0.8),
    log_density=log_density,
    burnin=200,
    proposal_scale=0.1,
    rng=Random.MersenneTwister(123)
)

println("MCMC results:")
println("  Parameter space: $(result.space.n_params) q-parameters")
println("  Block keys: ", result.space.block_keys)
println("  Acceptance rate: ", round(result.acceptance_rate; digits=3))
println("  Chain length: ", length(result.matrices))

# Posterior summary
posterior_ais = [assortativity_index(cm, :daytype) for cm in result.matrices]
mean_ai = sum(posterior_ais) / length(posterior_ais)
ai_std = sqrt(sum((ai - mean_ai)^2 for ai in posterior_ais) / length(posterior_ais))
println("  Mean posterior AI: ", round(mean_ai; digits=4))
println("  Target AI: ", target_ai)
println("  AI std: ", round(ai_std; digits=4))
```

## Step 8: Per-Block q-Parameters

Different age blocks may have different mixing preferences. We can specify
per-block parameters to capture this heterogeneity.

```@example v06
# Per-block: children mix more assortatively by daytype than adults
block_params = Dict(
    (1,1) => BlockAssortativityParams(q=Dict(:daytype => 0.6)),  # children x children
    (2,2) => BlockAssortativityParams(q=Dict(:daytype => 0.2)),  # adults x adults
    (3,3) => BlockAssortativityParams(q=Dict(:daytype => 0.1)),  # elderly x elderly
)
pspec_hetero = ParameterizedConstrainedLift(spec;
    block_params=block_params,
    default_params=BlockAssortativityParams(q=Dict(:daytype => 0.0))
)

cm_hetero = cm_age_sym ⊠ pspec_hetero
println("Heterogeneous q reconstruction:")
println("  Overall AI: ", round(assortativity_index(cm_hetero, :daytype); digits=4))
println("  rho: ", round(ρ(cm_hetero); digits=3))
```

## Summary

This vignette demonstrated the full categorical pipeline for partial-data
reconstruction using real POLYMOD data:

| Step | Operation | Categorical Interpretation |
|------|-----------|---------------------------|
| 1 | Source-stratified matrix | Section of survey functor (partial image) |
| 2 | Base matrix + symmetrise | Reciprocal object in age category |
| 3 | Constrained lift (`⊠`) | Section of coarsening morphism |
| 4 | Coarsen back (`↓`) | Verify lift is a section (left inverse) |
| 5 | q-parameters | Parameterized family of sections |
| 6 | Rejection sampling | Uniform exploration of fiber |
| 7 | MCMC | Targeted exploration with density on fiber |
| 8 | Per-block q | Heterogeneous sections |

The key insight: reconstruction from partial data is **not unique** -- it defines
a *fiber* over the coarsening morphism. The q-parameter family and MCMC sampler
let us explore this fiber under constraints (reciprocity, non-negativity) and
prior beliefs about mixing patterns.
