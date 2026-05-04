# ContACT.jl

*Applied category theory for structured contact matrices.*

## Overview

ContACT.jl provides a category-theoretic framework for constructing,
manipulating, and composing contact matrices over arbitrary finite survey
partitions: age bands, sex, region, occupation, risk groups, or products such as
age × sex. Built on [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl), it
provides formal guarantees (functoriality, idempotence, associativity) verified
both in tests and companion Lean 4 proofs.

## Algebraic Operators

All core operations have Unicode operators (type LaTeX name + TAB in the REPL):

| Operator | Input | Category-theoretic role |
|----------|-------|------------------------|
| `⊕` | `\oplus` | Commutative monoid (additive composition) |
| `⊗` | `\otimes` | Kronecker stratification functor |
| `↓` | `\downarrow` | Left Kan extension (coarsening) |
| `↑` | `\uparrow` | Parameterised refinement (with prior) |
| `⤊` | `\Uuparrow` | Activity refinement / hidden-stratum lift |
| `⊠` | `\boxtimes` | Generalized product lift |
| `▷` | `\triangleright` | Functor application (survey → matrix) |
| `∘` | `\circ` | Morphism composition (PartitionMap) |
| `↔` | `\leftrightarrow` | Reciprocity projection (symmetrisation) |
| `ρ` | `\rho` | Spectral radius (R₀ proxy) |

## Quick Example

```@example quickstart
using ContACT
using LinearAlgebra

# Define age groups and a contact matrix
partition = AgePartition([0, 18, 65])
pop = [1000.0, 3000.0, 500.0]
M = [2.0 1.0 0.5;
     1.0 3.0 1.0;
     0.5 1.0 1.5]
cm = ContactMatrix(M, partition, pop)
```

```@example quickstart
# Categorical partitions use symbols already present in survey tables
sex = CategoricalPartition(:sex;
    participant_col=:part_sex,
    contact_col=:cnt_sex,
    levels=["F", "M"],
    labels=["female", "male"],
)
group_labels(sex)
```

```@example quickstart
# Coarsen via ↓ (left Kan extension)
cm_coarse = cm ↓ AgePartition([0, 18])
matrix(cm_coarse)
```

```@example quickstart
# Compose maps with ∘, verify functoriality
fine = AgePartition([0, 18, 45, 65])
M4 = [2.0 1.0 0.5 0.2; 1.0 3.0 1.0 0.3; 0.5 1.0 2.5 0.8; 0.2 0.3 0.8 1.5]
pop4 = [1000.0, 2000.0, 1500.0, 500.0]
cm4 = ContactMatrix(M4, fine, pop4)

f = AgeMap(fine, partition)       # fine → 3 groups
g = AgeMap(partition, AgePartition([0, 18]))  # 3 → 2 groups
h = g ∘ f                        # composed: fine → 2 groups

# Functoriality: cm ↓ (g ∘ f) == (cm ↓ f) ↓ g
println("Functorial: $(matrix(cm4 ↓ h) ≈ matrix((cm4 ↓ f) ↓ g))")
```

```@example quickstart
# Symmetrise via ↔ (reciprocity projection)
cm_reciprocal = ↔(cm)
matrix(cm_reciprocal)
```

```@example quickstart
# Activity refinement uses respondent activity to lift a reciprocal matrix
# to partition × activity. Here `survey` would be a ContactSurvey over `partition`.
println("Activity lift operator: \\Uuparrow<TAB> gives ⤊")
```

```@example quickstart
# Generalized lifts construct age × SES/contact matrices from assumptions.
println("Generalized lift operator: \\boxtimes<TAB> gives ⊠")
```

```@example quickstart
# Spectral radius ρ — proportional to R₀
println("ρ(M) = $(round(ρ(cm); digits=2))")
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/epirecipes/ContACT.jl")
```

## Contents

```@contents
Pages = [
    "tutorials/getting_started.md",
    "tutorials/composition.md",
    "tutorials/categorical.md",
    "vignettes/04-generalized-contact-matrices/index.md",
    "api/types.md",
    "api/survey.md",
    "api/operations.md",
    "api/schemas.md",
]
```
