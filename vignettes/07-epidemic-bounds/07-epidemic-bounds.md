# Epidemic Bounds from Partial NGM Information
Simon Frost
2026-06-11

- [Overview](#overview)
- [Setup](#setup)
- [1. Scalar Final-Size Equation](#1-scalar-final-size-equation)
- [2. Multitype Final-Size Equation](#2-multitype-final-size-equation)
- [3. R₀ Bounds from Row/Column Sums](#3-r₀-bounds-from-rowcolumn-sums)
- [4. Detailed-Balance Bounds](#4-detailed-balance-bounds)
- [5. Final-Size Bounds](#5-final-size-bounds)
- [6. Total Final-Size Bounds](#6-total-final-size-bounds)
- [7. Using ContactMatrix Objects](#7-using-contactmatrix-objects)
- [8. Epidemic Uncertainty Across the q-Parameter
  Fiber](#8-epidemic-uncertainty-across-the-q-parameter-fiber)
- [9. Sensitivity to
  Transmissibility](#9-sensitivity-to-transmissibility)
- [10. Interpretation and Connection to
  ContACT](#10-interpretation-and-connection-to-contact)
- [Summary](#summary)

## Overview

This vignette demonstrates how to compute **bounds on R₀ and epidemic
final size** when the next-generation matrix (NGM) is only partially
known. The methods implement results from Britton, Poletti, Scarpa &
Pellis (2025, arXiv:2602.23885v2).

The key insight: when a contact matrix is known only up to its marginals
(row or column sums of the NGM), the spectral radius R₀ and the final
epidemic size are bounded — but not uniquely determined. ContACT.jl can
compute these **sharp** bounds analytically, and also explore the range
of epidemic outcomes across the **q-parameter fiber** of valid contact
matrices sharing the same marginals.

## Setup

``` julia
using ContACT
using LinearAlgebra
using Random
import ContACT: ×   # ContACT and LinearAlgebra both export ×; use ContACT's product
```

## 1. Scalar Final-Size Equation

The classic final-size relation for a homogeneous SIR model states:

$$1 - \tau = e^{-R_0 \cdot \tau}$$

where τ is the fraction of the population ultimately infected.

``` julia
# Subcritical: R₀ ≤ 1 → no epidemic
τ_sub = solve_final_size_scalar(0.8)
println("R₀ = 0.8 → τ = $τ_sub")

# Supercritical examples
for R0_val in [1.5, 2.0, 3.0, 5.0]
    τ = solve_final_size_scalar(R0_val)
    println("R₀ = $R0_val → τ = $(round(τ; digits=4))")
end
```

    R₀ = 0.8 → τ = 0.0
    R₀ = 1.5 → τ = 0.5828
    R₀ = 2.0 → τ = 0.7968
    R₀ = 3.0 → τ = 0.9405
    R₀ = 5.0 → τ = 0.993

## 2. Multitype Final-Size Equation

For a structured population with NGM `K` and population fractions `π`,
the final-size equations become coupled:

$$\tau_i = 1 - \exp\!\left(-\sum_j K_{ij}\,\pi_j\,\tau_j\right)$$

``` julia
# Two age groups with asymmetric mixing
K = [2.0 0.5;
     0.3 1.8]
π = [0.4, 0.6]

τ = solve_final_size_vector(K, π)
println("Per-type final sizes: τ₁ = $(round(τ[1]; digits=4)), τ₂ = $(round(τ[2]; digits=4))")
τ_bar = sum(π .* τ)
println("Population-average final size: τ̄ = $(round(τ_bar; digits=4))")
```

    Per-type final sizes: τ₁ = 0.8823, τ₂ = 0.8339
    Population-average final size: τ̄ = 0.8533

## 3. R₀ Bounds from Row/Column Sums

When we know only the **row sums** rᵢ = Σⱼ Kᵢⱼ (total reproductive
output of type i) or **column sums** cⱼ = Σᵢ Kᵢⱼ (total reproductive
input to type j), Theorem 3.1 gives sharp bounds:

- Row sums known: min(rᵢ) ≤ R₀ ≤ max(rᵢ)
- Column sums known: min(cⱼ) ≤ R₀ ≤ max(cⱼ)

``` julia
# A 4-type population (e.g., children, young adults, adults, elderly)
K4 = [3.0 1.0 0.5 0.1;
      0.8 2.5 1.2 0.3;
      0.4 1.0 2.0 0.6;
      0.1 0.3 0.5 1.0]
π4 = [0.2, 0.3, 0.35, 0.15]

# Actual R₀
R0_actual = maximum(abs.(eigvals(K4)))
println("Actual R₀ = $(round(R0_actual; digits=4))")

# Bounds from row sums
b_row = r0_bounds(K4; info=:row)
println("Row-sum bounds: $b_row")

# Bounds from column sums
b_col = r0_bounds(K4; info=:col)
println("Col-sum bounds: $b_col")

# Bounds from both
b_both = r0_bounds(K4; info=:both)
println("Combined bounds: $b_both")
```

    Actual R₀ = 4.2665
    Row-sum bounds: [1.9, 4.8]
    Col-sum bounds: [2.0, 4.8]
    Combined bounds: [2.0, 4.8]

## 4. Detailed-Balance Bounds

When the contact matrix satisfies **reciprocity** (πᵢKᵢⱼ = πⱼKⱼᵢ),
Theorem 3.2 provides tighter bounds using the weighted Cauchy-Schwarz
inequality:

``` julia
# Symmetric contact matrix (reciprocal by construction with equal populations)
K_sym = [3.0 1.2 0.4;
         1.2 2.5 0.8;
         0.4 0.8 1.5]
π_sym = [1/3, 1/3, 1/3]

R0_sym = maximum(abs.(eigvals(K_sym)))
println("Actual R₀ = $(round(R0_sym; digits=4))")

b_general = r0_bounds(K_sym; info=:row)
b_db = r0_bounds_detailed_balance(K_sym, π_sym; info=:row)
println("General bounds:          $b_general")
println("Detailed-balance bounds: $b_db")
println("Improvement: interval width $(round(b_general.upper - b_general.lower; digits=4)) → $(round(b_db.upper - b_db.lower; digits=4))")
```

    Actual R₀ = 4.2242
    General bounds:          [2.7, 4.6000000000000005]
    Detailed-balance bounds: [4.029061098237818, 4.6000000000000005]
    Improvement: interval width 1.9 → 0.5709

## 5. Final-Size Bounds

Beyond R₀, Theorems 3.3-3.4 bound the **per-type final sizes** τᵢ using
knowledge of row or column sums:

``` julia
# Using the 4-type matrix
τ_actual = solve_final_size_vector(K4, π4)
println("Actual per-type final sizes:")
for (i, τi) in enumerate(τ_actual)
    println("  Type $i: $(round(τi; digits=4))")
end

# Bounds from column sums (Theorem 3.3)
fs_col = final_size_bounds(K4, π4; info=:col)
println("\nColumn-sum bounds on τᵢ:")
for i in 1:4
    println("  Type $i: [$(round(fs_col.lower[i]; digits=4)), $(round(fs_col.upper[i]; digits=4))] (actual: $(round(τ_actual[i]; digits=4)))")
end

# Bounds from row sums (Theorem 3.4)
fs_row = final_size_bounds(K4, π4; info=:row)
println("\nRow-sum bounds on τᵢ:")
for i in 1:4
    println("  Type $i: [$(round(fs_row.lower[i]; digits=4)), $(round(fs_row.upper[i]; digits=4))] (actual: $(round(τ_actual[i]; digits=4)))")
end
```

    Actual per-type final sizes:
      Type 1: 0.9926
      Type 2: 0.9878
      Type 3: 0.9679
      Type 4: 0.9518

    Column-sum bounds on τᵢ:
      Type 1: [0.9234, 0.9994] (actual: 0.9926)
      Type 2: [0.8523, 0.996] (actual: 0.9878)
      Type 3: [0.7617, 0.984] (actual: 0.9679)
      Type 4: [0.7968, 0.9899] (actual: 0.9518)

    Row-sum bounds on τᵢ:
      Type 1: [0.0, 1.0] (actual: 0.9926)
      Type 2: [0.0, 0.9997] (actual: 0.9878)
      Type 3: [0.0, 0.9993] (actual: 0.9679)
      Type 4: [0.0, 1.0] (actual: 0.9518)

## 6. Total Final-Size Bounds

Theorem 3.5 bounds the **population-average** final size τ̄ = Σᵢ πᵢτᵢ:

``` julia
τ_bar_actual = sum(π4 .* τ_actual)
println("Actual total final size: τ̄ = $(round(τ_bar_actual; digits=4))")

tfs_row = total_final_size_bounds(K4, π4; info=:row)
println("Row-sum bounds on τ̄: $tfs_row")

tfs_col = total_final_size_bounds(K4, π4; info=:col)
println("Col-sum bounds on τ̄: $tfs_col")
```

    Actual total final size: τ̄ = 0.9764
    Row-sum bounds on τ̄: [0.11508653837950872, 0.9811025384344483]
    Col-sum bounds on τ̄: [0.8264886099666194, 0.9915302492004117]

## 7. Using ContactMatrix Objects

All epidemic bounds functions have convenience wrappers for
`ContactMatrix` inputs. These automatically compute the NGM and handle
the transpose convention:

``` julia
# Build a realistic 3-age-group contact matrix (POLYMOD-like).
# Symmetrise so total contacts are reciprocal — the detailed-balance bound
# (Theorem 3.2) is only valid when the NGM satisfies detailed balance.
part = AgePartition([0, 18, 65])
pop = [11_000.0, 33_000.0, 9_500.0]
M = [7.0 2.5 1.0;
     2.0 8.0 2.0;
     0.5 2.0 4.0]
cm = ↔(ContactMatrix(M, part, pop))

println("Contact matrix (per-capita rates):")
display(matrix(cm))

# R₀ bounds
R0_cm = R₀(cm)
println("\nActual R₀ = $(round(R0_cm; digits=4))")
b = r0_bounds(cm; info=:row)
println("Row-sum bounds: $b")
b_db = r0_bounds_detailed_balance(cm; info=:row)
println("Detailed-balance bounds: $b_db")

# Final size bounds
fs = final_size_bounds(cm; info=:col)
println("\nFinal size bounds (col sums):")
for (i, lbl) in enumerate(group_labels(cm))
    println("  $lbl: [$(round(fs.lower[i]; digits=4)), $(round(fs.upper[i]; digits=4))]")
end

# Total final size
tfs = total_final_size_bounds(cm; info=:row)
println("\nTotal final size bounds: $tfs")
```

    Contact matrix (per-capita rates):

    Actual R₀ = 11.0617
    Row-sum bounds: [9.263157894736842, 12.431818181818182]
    Detailed-balance bounds: [10.950403380600031, 12.431818181818182]

    Final size bounds (col sums):
      [0,18): [0.9997, 1.0]
      [18,65): [0.9929, 1.0]
      65+: [0.9974, 1.0]

    Total final size bounds: [0.1775532351727194, 0.9999816584579502]

    3×3 Matrix{Float64}:
     7.0       1.58333  0.789474
     4.75      8.0      4.47368
     0.681818  1.28788  4.0

## 8. Epidemic Uncertainty Across the q-Parameter Fiber

The most powerful application combines epidemic bounds with ContACT’s
**constrained lift** machinery. Different contact matrices sharing the
same marginals (the “fiber” of valid reconstructions) yield different
epidemic outcomes:

``` julia
Random.seed!(42)

# Start with a coarse age-only matrix
age_part = AgePartition([0, 30, 60])
age_pop = [15_000.0, 20_000.0, 10_000.0]
M_age = [5.0 2.0 0.5;
         2.0 6.0 1.5;
         0.5 1.5 3.0]
cm_age = ↔(ContactMatrix(M_age, age_part, age_pop))

# Refine by activity level (high/low) - the 'true' distribution is unknown
activity = CategoricalPartition(:activity, ["high", "low"])
prod_part = age_part × activity
prod_pop = [6000.0, 9000.0, 8000.0, 12000.0, 3000.0, 7000.0]

# Build an intermediate source-stratified matrix (source = age × activity, target = age)
# representing observed participant-side data before full reconstruction
n_age = n_groups(age_part)
n_prod = n_groups(prod_part)
M_inter = zeros(n_age, n_prod)
for j in 1:n_prod
    age_idx = ((j - 1) ÷ 2) + 1            # age group this source column belongs to
    for i in 1:n_age
        M_inter[i, j] = matrix(cm_age)[i, age_idx] *
            prod_pop[j] / sum(prod_pop[2*(age_idx-1)+1:2*age_idx])
    end
end
intermediate = SourceStratifiedContactMatrix(M_inter, age_part, prod_part, prod_pop)

# Set up the constrained lift and the base (proportional) reconstruction
spec = ConstrainedGeneralizedLift(intermediate)
full_cm = cm_age ⊠ spec

# Sample from the q-parameter fiber
samples = sample_constrained_lifts(cm_age, spec, 100; bounds=(-0.5, 0.5))
matrices = [s[2] for s in samples]

# Epidemic uncertainty across the fiber
epi = epidemic_uncertainty(matrices)
println("R₀ range across fiber: $(round(epi.r0.lower; digits=4)) – $(round(epi.r0.upper; digits=4))")
println("Final size range:      $(round(epi.final_size.lower; digits=4)) – $(round(epi.final_size.upper; digits=4))")

# Compare with analytic bounds (which don't require sampling)
b_analytic = r0_bounds(full_cm; info=:row)
println("\nAnalytic R₀ bounds (row sums): $b_analytic")
println("Fiber exploration is tighter: $(round(epi.r0.upper - epi.r0.lower; digits=4)) vs $(round(b_analytic.upper - b_analytic.lower; digits=4))")
```

    R₀ range across fiber: 8.3538 – 8.5207
    Final size range:      0.996 – 0.9961

    Analytic R₀ bounds (row sums): [3.038793103448276, 10.240384615384615]
    Fiber exploration is tighter: 0.1669 vs 7.2016

## 9. Sensitivity to Transmissibility

The bounds can be computed for different transmissibility parameters
(β/γ), showing how uncertainty in within-group mixing translates to
uncertainty in epidemic outcomes at different transmission intensities:

``` julia
println("Transmissibility │ R₀ bounds (row) │ Total final size bounds")
println("─────────────────┼─────────────────┼────────────────────────")
for β in [0.5, 1.0, 1.5, 2.0]
    b = r0_bounds(cm_age; info=:row, transmissibility=β)
    tfs = total_final_size_bounds(cm_age; info=:row, transmissibility=β)
    println("  β = $(rpad(β, 15)) │ [$(round(b.lower; digits=2)), $(round(b.upper; digits=2))]$(repeat(" ", 8 - length(string(round(b.upper; digits=2))))) │ [$(round(tfs.lower; digits=4)), $(round(tfs.upper; digits=4))]")
end
```

    Transmissibility │ R₀ bounds (row) │ Total final size bounds
    ─────────────────┼─────────────────┼────────────────────────
      β = 0.5             │ [2.94, 4.44]     │ [0.208, 0.9783]
      β = 1.0             │ [5.88, 8.88]     │ [0.2216, 0.9996]
      β = 1.5             │ [8.81, 13.31]    │ [0.2222, 1.0]
      β = 2.0             │ [11.75, 17.75]    │ [0.2222, 1.0]

## 10. Interpretation and Connection to ContACT

The epidemic bounds framework connects naturally to ContACT.jl’s core
ideas:

1.  **Coarsening loses information**: When a fine-grained contact matrix
    is coarsened (via `↓`), the within-block structure is lost. The
    epidemic bounds quantify how much this information loss affects R₀
    and final size predictions.

2.  **The q-parameter fiber**: All matrices sharing the same marginals
    form a fiber. The `epidemic_uncertainty` function explores this
    fiber to find the actual range of epidemic outcomes, which is
    generally **tighter** than the analytic bounds (which consider *all
    possible* matrices, not just those reachable through the
    q-parameterization).

3.  **Practical epidemiology**: When survey data provides only aggregate
    contact rates by age group, the within-group heterogeneity (by
    activity level, socioeconomic status, etc.) is unknown. These bounds
    tell us how much this ignorance matters for epidemic predictions.

## Summary

| Function                     | What it bounds  | Input needed                |
|------------------------------|-----------------|-----------------------------|
| `r0_bounds`                  | R₀              | Row/column sums of NGM      |
| `r0_bounds_detailed_balance` | R₀ (tighter)    | Row/col sums + reciprocity  |
| `final_size_bounds`          | Per-type τᵢ     | Row or column sums + π      |
| `total_final_size_bounds`    | Population τ̄    | Row or column sums + π      |
| `epidemic_uncertainty`       | R₀ and τ̄ ranges | Sampled matrices from fiber |
