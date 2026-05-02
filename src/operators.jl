"""
Unicode operators for ContACT.jl.

These provide a concise algebraic syntax for contact matrix operations.
Type the LaTeX name followed by TAB in the Julia REPL to enter these:
  \\oplus<TAB>  → ⊕
  \\otimes<TAB> → ⊗
  \\downarrow<TAB> → ↓
"""

"""
    a ⊕ b

Additive composition of contact matrices. Type `\\oplus<TAB>`.

Equivalent to `compose_matrices(a, b)`.

# Example
```julia
total = cm_home ⊕ cm_work ⊕ cm_school ⊕ cm_other
```
"""
⊕(a::ContactMatrix, b::ContactMatrix) = compose_matrices(a, b)

"""
    cm ⊗ coupling

Stratify a contact matrix by a coupling matrix. Type `\\otimes<TAB>`.

Equivalent to `stratify(cm, coupling)`.

# Example
```julia
# 3 regions with uniform mixing
spatial = cm ⊗ ones(3, 3) / 3
```
"""
⊗(cm::ContactMatrix, coupling::AbstractMatrix) = stratify(cm, coupling)

"""
    cm ↓ coarse_partition

Coarsen a contact matrix to a coarser age partition. Type `\\downarrow<TAB>`.

Equivalent to `coarsen(cm, coarse_partition)`.

# Example
```julia
coarse = AgePartition([0, 18, 65])
cm_coarse = cm ↓ coarse
```
"""
↓(cm::ContactMatrix, coarse::AgePartition) = coarsen(cm, coarse)
