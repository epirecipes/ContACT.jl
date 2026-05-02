# ContACT.jl

*Applied category theory for contact matrices in infectious disease modelling.*

## Overview

ContACT.jl provides a category-theoretic framework for constructing,
manipulating, and composing age-structured contact matrices from social mixing
survey data. Built on [Catlab.jl](https://github.com/AlgebraicJulia/Catlab.jl),
it provides formal guarantees (functoriality, idempotence, associativity) verified
both in tests and companion Lean 4 proofs.

## Key Features

| Operation | Operator | Category-theoretic interpretation |
|-----------|----------|----------------------------------|
| Composition | `⊕` | Commutative monoid on Hom(G,G) |
| Stratification | `⊗` | Kronecker product with coupling |
| Coarsening | `↓` | Left Kan extension along surjection |
| Symmetrisation | `symmetrise` | Idempotent endomorphism |
| Survey → Matrix | `compute_matrix` | Restricted functor |

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
# Symmetrise (idempotent endomorphism)
cm_sym = symmetrise(cm)
matrix(cm_sym)
```

```@example quickstart
# Coarsen to 2 groups (left Kan extension)
coarse = AgePartition([0, 18])
cm_coarse = cm ↓ coarse
matrix(cm_coarse)
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/your-org/ContACT.jl")
```

## Contents

```@contents
Pages = [
    "tutorials/getting_started.md",
    "tutorials/composition.md",
    "tutorials/categorical.md",
    "api/types.md",
    "api/survey.md",
    "api/operations.md",
    "api/schemas.md",
]
```
