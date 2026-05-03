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
| `▷` | `\triangleright<TAB>` | Functor application |
| `∘` | `\circ<TAB>` | Partition map composition |
| `↔` | `\leftrightarrow<TAB>` | Symmetrisation / reciprocity |
| `ρ` | `\rho<TAB>` | Spectral radius |

```@docs
⊕
⊗
↓
↑
▷
↔
ρ
RefinementPrior
```

## Map Composition

```@docs
Base.:∘(::PartitionMap{D, P, Q} where {P<:AbstractPartition{D}, Q<:AbstractPartition{D}}, ::PartitionMap{D, P, Q} where {P<:AbstractPartition{D}, Q<:AbstractPartition{D}}) where D
```

## Coarsening & Refinement

```@docs
PartitionMap
AgeMap
coarsen
refine
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
