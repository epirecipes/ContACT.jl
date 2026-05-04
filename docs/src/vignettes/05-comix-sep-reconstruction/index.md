# CoMix/SEP Partial-Data Reconstruction

This vignette demonstrates how to reconstruct a fully stratified contact
matrix from partial survey data, following the approach of Di Domenico,
Reichmuth & Althaus (2025). The key insight is that participant-side
covariates (age, socioeconomic position, education) are fully observed, but
contact-side covariates beyond age are missing.

## Categorical Framework

The reconstruction can be understood categorically:

1. **Observed data** lives in a *source-stratified* space: contacts are
   classified only by target age, but participants are classified by the full
   product partition `age × SEP × education`.
2. **The base matrix** over age provides a reciprocal constraint: total
   contacts between age groups must be balanced.
3. **Reconstruction** is an assumption-driven *section* of the coarsening
   morphism `age × SEP × edu → age`. It is not a unique inverse — a family of
   valid reconstructions exists, parameterized by within-block assortativity.
4. **Feasible sampling** explores this family by varying q-parameters and
   rejecting non-positive solutions.

## Setup

```julia
using ContACT
import ContACT: ×

# Define partitions
age = IntervalPartition{:age,Float64}([0.0, 18.0, 65.0])  # 3 age groups
sep = CategoricalPartition(:sep, ["low_SEP", "high_SEP"])
edu = CategoricalPartition(:edu, ["low_edu", "high_edu"])

# Full product partition: 3 × 2 × 2 = 12 groups
full = age × sep × edu
```

## Source-Stratified Intermediate Matrix

In practice, this matrix is computed from survey data via
`compute_source_stratified_matrix`. Here we use a synthetic example:

```julia
# 3 target (age) rows × 12 source (age×sep×edu) columns
# Population of each source group
pop12 = [80.0, 40.0, 60.0, 20.0,   # children: low_SEP×low_edu, low_SEP×high_edu, high_SEP×low_edu, high_SEP×high_edu
         100.0, 80.0, 90.0, 70.0,  # adults: ...
         50.0, 30.0, 40.0, 30.0]   # elderly: ...

# Mean contacts reported by each source group with each target age group
M_inter = [
    5.0 4.5 4.0 3.5  1.5 1.2 1.8 1.4  0.8 0.6 0.7 0.5;  # contacts in children
    2.0 1.8 2.2 1.9  4.0 3.8 3.5 3.2  1.5 1.2 1.3 1.0;  # contacts in adults
    0.5 0.4 0.6 0.5  1.0 0.9 1.2 1.1  3.0 2.8 2.5 2.2   # contacts in elderly
]

intermediate = SourceStratifiedContactMatrix(M_inter, age, full, pop12)
```

## Base Reciprocal Age Matrix

The base matrix over age must be reciprocal in total-contact space
(`C = M × diag(N)` is symmetric):

```julia
age_pop = [200.0, 340.0, 150.0]  # sum of source pops per age group
M_age = [4.5 2.0 0.6;
         1.18 3.6 1.2;
         0.8 2.72 2.7]
# Verify reciprocity: C = M * diag(N) should be approximately symmetric
base = ContactMatrix(M_age, age, age_pop)
```

## Constrained Lift Specification

```julia
source_to_age = PartitionMap(full, age)
spec = ConstrainedGeneralizedLift(intermediate; source_map=source_to_age)
```

## Proportionate Reconstruction (q = 0)

The baseline reconstruction assumes proportionate mixing within each age block:

```julia
full_cm = base ⊠ spec

# Verify invariants
@assert matrix(full_cm ↓ age) ≈ matrix(base)  # coarsens back to base
C = matrix(full_cm) * Diagonal(population(full_cm))
@assert C ≈ C'  # reciprocal
```

## Parameterized Reconstruction (q ≠ 0)

Adding assortativity parameters shifts mass toward same-group contacts:

```julia
# Global SEP assortativity
params = BlockAssortativityParams(q=Dict(:sep => 0.3, :edu => 0.1))
pspec = ParameterizedConstrainedLift(spec; default_params=params)
full_assort = base ⊠ pspec

# Assortativity index increases with q
ai_prop = assortativity_index(full_cm, :sep)
ai_assort = assortativity_index(full_assort, :sep)
@assert ai_assort > ai_prop
```

## Per-Block Parameters

Different age-block interactions can have different assortativity:

```julia
block_params = Dict(
    (1, 1) => BlockAssortativityParams(q=Dict(:sep => 0.5)),   # children mixing
    (2, 2) => BlockAssortativityParams(q=Dict(:sep => 0.2)),   # adult mixing
    (3, 3) => BlockAssortativityParams(q=Dict(:sep => 0.1)),   # elderly mixing
)
pspec_blocks = ParameterizedConstrainedLift(spec;
    block_params=block_params,
    default_params=BlockAssortativityParams())
full_blocks = base ⊠ pspec_blocks
```

## Feasible Sampling

Explore the space of valid reconstructions:

```julia
using Random
samples = sample_constrained_lifts(base, spec, 50;
    dimensions=[:sep, :edu],
    bounds=(-0.5, 0.5),
    rng=MersenneTwister(123))

# Analyze the sample
for (params, cm) in samples
    R0_val = basic_reproduction_number(cm)
    ai_sep = assortativity_index(cm, :sep)
    println("q_sep=$(params.q[:sep]) q_edu=$(params.q[:edu]) R0=$R0_val AI_sep=$ai_sep")
end
```

## Targeted Control Analysis

Compute type-reproduction numbers for targeted interventions:

```julia
# Target high-SEP groups for vaccination
# In the 12-group product partition, high-SEP groups are those with
# factor index 2 in the SEP dimension
Tg = type_reproduction_number(full_assort, [3, 4, 7, 8, 11, 12])
threshold = control_threshold(Tg)
effort = control_effort(full_assort, [3, 4, 7, 8, 11, 12], threshold)
println("Type reproduction number: $Tg")
println("Control threshold: $threshold")
println("Control effort: $effort")
```

## Marginal Matrices

Extract dimension-specific marginals from the full reconstruction:

```julia
# SEP-only marginal (2×2)
cm_sep = marginal_matrix(full_assort, :sep)

# Education-only marginal (2×2)
cm_edu = marginal_matrix(full_assort, :edu)

# Age-only (recovers base)
cm_age = marginal_matrix(full_assort, :age)
@assert matrix(cm_age) ≈ matrix(base)
```

## Summary

The CoMix/SEP reconstruction workflow in ContACT.jl:

| Step | Function | Operator |
|------|----------|----------|
| Build intermediate | `compute_source_stratified_matrix` | — |
| Specify constraints | `ConstrainedGeneralizedLift` | — |
| Proportionate lift | `constrained_generalize` | `⊠` |
| Parameterized lift | `ParameterizedConstrainedLift` | `⊠` |
| Feasibility check | `is_feasible` | — |
| Sample family | `sample_constrained_lifts` | — |
| SEP marginal | `marginal_matrix(cm, :sep)` | — |
| Assortativity | `assortativity_index(cm, :sep)` | — |
| Type repro. number | `type_reproduction_number` | — |
| Control threshold | `control_threshold` | — |
| Control effort | `control_effort` | — |
