# Matrix Operations

## Unicode Operators

ContACT.jl provides a concise algebraic syntax via Unicode operators.
Type the LaTeX name followed by TAB in the Julia REPL:

| Operator | Input | Operation |
|----------|-------|-----------|
| `⊕` | `\oplus<TAB>` | Additive composition |
| `⊗` | `\otimes<TAB>` | Stratification (Kronecker) |
| `↓` | `\downarrow<TAB>` | Coarsening (left Kan extension) |
| `↑` | `\uparrow<TAB>` | Refinement (with prior) |
| `⤊` | `\Uuparrow<TAB>` | Activity refinement / hidden-stratum lift |
| `⊠` | `\boxtimes<TAB>` | Generalized product lift |
| `▷` | `\triangleright<TAB>` | Functor application |
| `∘` | `\circ<TAB>` | Partition map composition |
| `↔` | `\leftrightarrow<TAB>` | Symmetrisation / reciprocity |
| `ρ` | `\rho<TAB>` | Spectral radius |

```@docs
⊕
⊗
↓
↑
⤊
⊠
▷
↔
ρ
RefinementPrior
```

## Map Composition

```@docs
Base.:∘(::PartitionMap, ::PartitionMap)
```

## Coarsening & Refinement

```@docs
PartitionMap
AgeMap
coarsen
refine
```

## Activity Refinement

```@docs
ActivityRefinement
ActivityMixingKernel
AssortativeMixing
DisassortativeMixing
ProportionateMixing
activity_partition
activity_mixing_plan
activity_refine
proportionate_mixing
```

## Generalized Contact Matrices

```@docs
GeneralizedLift
GeneralizedMixingKernel
RandomMixing
BlockMixing
AssortativeDimensionMixing
product_population
generalize
generalized_lift
```

## Partial-Data Reconstruction

```@docs
SourceStratifiedContactMatrix
target_partition
source_partition
n_target_groups
n_source_groups
target_group_labels
source_group_labels
compute_source_stratified_matrix
source_total_contacts
coarsen_sources
align_source_stratified_matrix
ConstrainedGeneralizedLift
full_partition
intermediate_matrix
source_map
structural_zeros
constrained_generalize
BlockAssortativityParams
ParameterizedConstrainedLift
is_feasible
sample_constrained_lifts
QParameterSpace
sample_perblock_lifts
mcmc_constrained_lifts
MCMCResult
```

## Composition

```@docs
compose_matrices
```

## Stratification

```@docs
stratify
```

## Symmetrisation

```@docs
symmetrise
```

## Reprojection

```@docs
reproject
```

## Population transport

```@docs
transport_population
```

## Unit semantics

`ContactCounts`, `MeanContacts` and `PerCapitaRate` are three representations of one
underlying object — total contacts `T[i,j] = M[i,j]·N_j`, rescaled by populations
along the contactee (row) and participant (column) indices. Conversions dispatch on
the semantics the matrix actually carries, and `symmetrise` is natural with respect
to them.

```@docs
population_exponents
reinterpret_units
to_mean_contacts
to_per_capita
to_counts
```

## Utilities

```@docs
spectral_radius
next_generation_matrix
basic_reproduction_number
calibrate_transmissibility
R0
R₀
marginal_matrix
assortativity_index
type_reproduction_number
control_threshold
control_effort
```

## Epidemic Bounds

Functions implementing R₀ and final-size bounds from partial NGM information
(Britton, Poletti, Scarpa & Pellis, 2025).

```@docs
EpidemicBounds
r0_bounds
r0_bounds_detailed_balance
final_size_bounds
total_final_size_bounds
solve_final_size_scalar
solve_final_size_ext
solve_final_size_vector
epidemic_uncertainty
```
