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

## Utilities

```@docs
to_per_capita
to_counts
spectral_radius
```
