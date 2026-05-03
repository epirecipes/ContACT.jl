"""
Unicode operators for ContACT.jl.

These provide a concise algebraic syntax for contact matrix operations.
Type the LaTeX name followed by TAB in the Julia REPL to enter these:

| Operator | LaTeX | Operation | Example |
|----------|-------|-----------|---------|
| `⊕` | `\\oplus` | Additive composition | `cm_home ⊕ cm_work` |
| `⊗` | `\\otimes` | Stratification (Kronecker) | `cm ⊗ coupling` |
| `↓` | `\\downarrow` | Coarsening | `cm ↓ coarse` |
| `↑` | `\\uparrow` | Refinement | `cm ↑ prior` |
| `⤊` | `\\Uuparrow` | Activity refinement | `cm ⤊ ActivityRefinement(survey)` |
| `▷` | `\\triangleright` | Functor (compute matrix) | `survey ▷ partition` |
| `∘` | `\\circ` | Map composition | `g ∘ f` |
| `↔` | `\\leftrightarrow` | Symmetrisation/reciprocity | `↔(cm)` |
| `ρ` | `\\rho` | Spectral radius | `ρ(cm)` |
"""

# ═══════════════════════════════════════════════════════════════════════════════
# Additive composition: ⊕
# ═══════════════════════════════════════════════════════════════════════════════

"""
    a ⊕ b

Additive composition of contact matrices (commutative monoid). Type `\\oplus<TAB>`.

Categorically, this is the monoidal product in the category of contact matrices
over a fixed age partition: (ContactMat, ⊕, 𝟎) where 𝟎 is the zero matrix.

# Properties (proven in Lean)
- Associativity: `(A ⊕ B) ⊕ C == A ⊕ (B ⊕ C)`
- Commutativity: `A ⊕ B == B ⊕ A`
- Identity: `A ⊕ 𝟎 == A`

# Example
```julia
cm_total = cm_home ⊕ cm_work ⊕ cm_school ⊕ cm_other
```
"""
⊕(a::ContactMatrix, b::ContactMatrix) = compose_matrices(a, b)

# ═══════════════════════════════════════════════════════════════════════════════
# Stratification (Kronecker product): ⊗
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cm ⊗ coupling

Stratify a contact matrix by a coupling matrix (Kronecker product). Type `\\otimes<TAB>`.

Creates a block-structured matrix where local age contacts are modulated by
inter-stratum coupling: `M_spatial = coupling ⊗ M_local`.

# Example
```julia
coupling = [0.8 0.2; 0.2 0.8]
cm_spatial = cm ⊗ coupling   # 2 regions × n ages
```
"""
⊗(cm::ContactMatrix, coupling::AbstractMatrix) = stratify(cm, coupling)

# ═══════════════════════════════════════════════════════════════════════════════
# Coarsening (left Kan extension): ↓
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cm ↓ coarse_partition
    cm ↓ f::PartitionMap

Coarsen a contact matrix along a surjective partition map. Type `\\downarrow<TAB>`.

This is the left Kan extension: it pushes forward contact data while preserving
total contacts. Satisfies functoriality (proven in Lean):

    `cm ↓ (g ∘ f) == (cm ↓ f) ↓ g`

# Example
```julia
cm_coarse = cm ↓ AgePartition([0, 18, 65])
cm_coarse = cm ↓ PartitionMap(fine, coarse)
```
"""
↓(cm::ContactMatrix, coarse::AbstractPartition) = coarsen(cm, coarse)
↓(cm::ContactMatrix, f::PartitionMap) = coarsen(cm, f)

# ═══════════════════════════════════════════════════════════════════════════════
# Refinement (parameterised disaggregation): ↑
# ═══════════════════════════════════════════════════════════════════════════════

"""
    RefinementPrior(partition, population)

Packages a target fine partition with the distributional prior (population per group)
needed for proportional refinement. Use with the `↑` operator.

# Example
```julia
prior = RefinementPrior(AgePartition([0, 5, 18, 65]), [800, 1200, 3000, 500])
cm_fine = cm ↑ prior
```
"""
struct RefinementPrior
    partition::AbstractPartition
    population::Vector{Float64}

    function RefinementPrior(partition::AbstractPartition, population::AbstractVector{<:Real})
        n_groups(partition) == length(population) || throw(DimensionMismatch(
            "partition has $(n_groups(partition)) groups but population has $(length(population)) entries"))
        all(x -> isfinite(x) && x >= 0, population) ||
            throw(ArgumentError("population entries must be finite and non-negative"))
        new(partition, Float64.(population))
    end
end

"""
    cm ↑ prior::RefinementPrior
    cm ↑ spec::ActivityRefinement

Refine a coarse contact matrix to a finer partition. Type `\\uparrow<TAB>`.

Unlike coarsening (which is canonical), refinement requires auxiliary assumptions
encoded in a `RefinementPrior`. This is NOT an inverse of `↓`.

When the right-hand side is an `ActivityRefinement`, this performs the
Britton-Ball style activity lift. Use `⤊` for a visually distinct alias.

# Example
```julia
prior = RefinementPrior(AgePartition([0, 5, 18, 65]), fine_pop)
cm_fine = cm ↑ prior
```
"""
↑(cm::ContactMatrix, prior::RefinementPrior) = refine(cm, prior.partition, prior.population)
↑(cm::ContactMatrix, spec::ActivityRefinement) = activity_refine(cm, spec)

"""
    cm ⤊ spec::ActivityRefinement

Refine a reciprocal contact matrix by respondent activity strata. Type
`\\Uuparrow<TAB>`.

The double upward arrow denotes a lift to the product partition
`cm.partition × activity`, with the activity coupling specified by `spec`.

# Example
```julia
spec = ActivityRefinement(survey; n=2, mixing=:proportionate)
cm_activity = ↔(cm) ⤊ spec
```
"""
⤊(cm::ContactMatrix, spec::ActivityRefinement) = activity_refine(cm, spec)

# ═══════════════════════════════════════════════════════════════════════════════
# Functor application: ▷
# ═══════════════════════════════════════════════════════════════════════════════

"""
    survey ▷ partition

Apply the survey-to-matrix functor. Type `\\triangleright<TAB>`.

Categorically, this is the image of `survey` under the functor
F: ContactSurvey → ContactMatrix. For the full API with population and
weights, use `compute_matrix(survey, partition; kwargs...)`.

# Example
```julia
cm = survey ▷ partition
```
"""
▷(survey::ContactSurvey, partition::AbstractPartition) = compute_matrix(survey, partition)

# ═══════════════════════════════════════════════════════════════════════════════
# PartitionMap composition: ∘
# ═══════════════════════════════════════════════════════════════════════════════

"""
    g ∘ f

Compose partition maps (right-to-left, standard mathematical convention).
Type `\\circ<TAB>`.

Given `f: A → B` and `g: B → C`, yields `g ∘ f: A → C`.
Satisfies functoriality (proven in Lean):

    `coarsen(cm, g ∘ f) == coarsen(coarsen(cm, f), g)`

# Example
```julia
f = AgeMap(fine, medium)      # fine → medium
g = AgeMap(medium, coarse)    # medium → coarse
h = g ∘ f                     # fine → coarse (composed)
```
"""
function Base.:∘(g::PartitionMap, f::PartitionMap)
    same_partition(f.codomain, g.domain) || throw(ArgumentError(
        "Cannot compose: f codomain does not match g domain"))
    # Compose the underlying FinFunctions
    f_assignments = collect(f.mapping)
    g_assignments = collect(g.mapping)
    composed = [g_assignments[f_assignments[i]] for i in eachindex(f_assignments)]
    PartitionMap(f.domain, g.codomain, composed)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Symmetrisation / reciprocity: ↔
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ↔(cm)

Symmetrise a contact matrix by enforcing reciprocal total contacts.
Type `\\leftrightarrow<TAB>`.

The arrow denotes bidirectional balance:

    `matrix(↔(cm))[i,j] * N[j] == matrix(↔(cm))[j,i] * N[i]`

This is a prefix operator alias for `symmetrise(cm)`.

# Example
```julia
cm_reciprocal = ↔(cm)
```
"""
↔(cm::ContactMatrix) = symmetrise(cm)

# ═══════════════════════════════════════════════════════════════════════════════
# Spectral radius: ρ
# ═══════════════════════════════════════════════════════════════════════════════

"""
    ρ(cm)

Spectral radius (dominant eigenvalue) of a contact matrix. Type `\\rho<TAB>`.

In epidemiology, `ρ(M)` is proportional to R₀ for frequency-dependent
transmission: R₀ = β/γ · ρ(M).

# Example
```julia
ρ(cm_total)           # basic reproductive ratio proxy
ρ(cm_lockdown) / ρ(cm_total)  # relative reduction
```
"""
ρ(cm::ContactMatrix) = spectral_radius(cm)
