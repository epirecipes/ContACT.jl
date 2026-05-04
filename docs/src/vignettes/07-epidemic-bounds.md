# Epidemic Bounds from Partial NGM Information

## Overview

This vignette demonstrates how to compute **bounds on R₀ and epidemic final
size** when the next-generation matrix (NGM) is only partially known. The
methods implement results from Britton, Poletti, Scarpa & Pellis (2025,
arXiv:2602.23885v2).

The key insight: when a contact matrix is known only up to its marginals (row
or column sums of the NGM), the spectral radius R₀ and the final epidemic size
are bounded — but not uniquely determined. ContACT.jl can compute these
**sharp** bounds analytically, and also explore the range of epidemic outcomes
across the **q-parameter fiber** of valid contact matrices sharing the same
marginals.

## Setup

```@example v07
using ContACT
using LinearAlgebra: eigvals, Diagonal
import ContACT: ×
using Random
```

## 1. Scalar Final-Size Equation

The classic final-size relation for a homogeneous SIR model states:

```math
1 - \tau = e^{-R_0 \cdot \tau}
```

where τ is the fraction of the population ultimately infected.

```@example v07
# Subcritical: R₀ ≤ 1 → no epidemic
τ_sub = solve_final_size_scalar(0.8)
println("R₀ = 0.8 → τ = $τ_sub")

# Supercritical examples
for R0_val in [1.5, 2.0, 3.0, 5.0]
    τ = solve_final_size_scalar(R0_val)
    println("R₀ = $R0_val → τ = $(round(τ; digits=4))")
end
```

## 2. Multitype Final-Size Equation

For a structured population with NGM `K` and population fractions `π`,
the final-size equations become coupled:

```math
\tau_i = 1 - \exp\!\left(-\sum_j K_{ij}\,\pi_j\,\tau_j\right)
```

```@example v07
# Two age groups with asymmetric mixing
K = [2.0 0.5;
     0.3 1.8]
π = [0.4, 0.6]

τ = solve_final_size_vector(K, π)
println("Per-type final sizes: τ₁ = $(round(τ[1]; digits=4)), τ₂ = $(round(τ[2]; digits=4))")
τ_bar = sum(π .* τ)
println("Population-average final size: τ̄ = $(round(τ_bar; digits=4))")
```

## 3. R₀ Bounds from Row/Column Sums

When we know only the **row sums** rᵢ = Σⱼ Kᵢⱼ (total reproductive output of
type i) or **column sums** cⱼ = Σᵢ Kᵢⱼ (total reproductive input to type j),
Theorem 3.1 gives sharp bounds:

- Row sums known: min(rᵢ) ≤ R₀ ≤ max(rᵢ)
- Column sums known: min(cⱼ) ≤ R₀ ≤ max(cⱼ)

```@example v07
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

## 4. Detailed-Balance Bounds

When the contact matrix satisfies **reciprocity** (πᵢKᵢⱼ = πⱼKⱼᵢ), Theorem 3.2
provides tighter bounds using the weighted Cauchy-Schwarz inequality:

```@example v07
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

## 5. Final-Size Bounds

Beyond R₀, Theorems 3.3-3.4 bound the **per-type final sizes** τᵢ using
knowledge of row or column sums:

```@example v07
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

## 6. Total Final-Size Bounds

Theorem 3.5 bounds the **population-average** final size τ̄ = Σᵢ πᵢτᵢ:

```@example v07
τ_bar_actual = sum(π4 .* τ_actual)
println("Actual total final size: τ̄ = $(round(τ_bar_actual; digits=4))")

tfs_row = total_final_size_bounds(K4, π4; info=:row)
println("Row-sum bounds on τ̄: $tfs_row")

tfs_col = total_final_size_bounds(K4, π4; info=:col)
println("Col-sum bounds on τ̄: $tfs_col")
```

## 7. Using ContactMatrix Objects

All epidemic bounds functions have convenience wrappers for `ContactMatrix`
inputs. These automatically compute the NGM and handle the transpose convention:

```@example v07
# Build a realistic 3-age-group contact matrix (POLYMOD-like)
part = AgePartition([0, 18, 65])
pop = [11_000.0, 33_000.0, 9_500.0]
M = [7.0 2.5 1.0;
     2.0 8.0 2.0;
     0.5 2.0 4.0]
cm = ContactMatrix(M, part, pop)

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

## 8. Epidemic Uncertainty Across the q-Parameter Fiber

The most powerful application combines epidemic bounds with ContACT's
**constrained lift** machinery. Different contact matrices sharing the same
marginals (the "fiber" of valid reconstructions) yield different epidemic
outcomes:

```@example v07
Random.seed!(42)

# Start with a coarse age-only matrix (must be reciprocal)
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

# Build an intermediate source-stratified matrix (source=age×activity, target=age)
# This represents observed participant-side data before full reconstruction
n_age = n_groups(age_part)
n_prod = n_groups(prod_part)
M_inter = zeros(n_age, n_prod)
for j in 1:n_prod
    age_idx = ((j - 1) ÷ 2) + 1  # which age group this source belongs to
    for i in 1:n_age
        M_inter[i, j] = matrix(cm_age)[i, age_idx] * prod_pop[j] / sum(prod_pop[2*(age_idx-1)+1:2*age_idx])
    end
end
intermediate = SourceStratifiedContactMatrix(M_inter, age_part, prod_part, prod_pop)

# Set up the constrained lift
spec = ConstrainedGeneralizedLift(intermediate)

# Compute the base (proportional) lift
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

## 9. Sensitivity to Transmissibility

The bounds can be computed for different transmissibility parameters (β/γ),
showing how uncertainty in within-group mixing translates to uncertainty in
epidemic outcomes at different transmission intensities:

```@example v07
println("Transmissibility │ R₀ bounds (row) │ Total final size bounds")
println("─────────────────┼─────────────────┼────────────────────────")
for β in [0.5, 1.0, 1.5, 2.0]
    b = r0_bounds(cm_age; info=:row, transmissibility=β)
    tfs = total_final_size_bounds(cm_age; info=:row, transmissibility=β)
    println("  β = $(rpad(β, 15)) │ [$(round(b.lower; digits=2)), $(round(b.upper; digits=2))]$(repeat(" ", 8 - length(string(round(b.upper; digits=2))))) │ [$(round(tfs.lower; digits=4)), $(round(tfs.upper; digits=4))]")
end
```

## 10. Interpretation and Connection to ContACT

The epidemic bounds framework connects naturally to ContACT.jl's core ideas:

1. **Coarsening loses information**: When a fine-grained contact matrix is
   coarsened (via `↓`), the within-block structure is lost. The epidemic bounds
   quantify how much this information loss affects R₀ and final size predictions.

2. **The q-parameter fiber**: All matrices sharing the same marginals form a
   fiber. The `epidemic_uncertainty` function explores this fiber to find the
   actual range of epidemic outcomes, which is generally **tighter** than the
   analytic bounds (which consider *all possible* matrices, not just those
   reachable through the q-parameterization).

3. **Practical epidemiology**: When survey data provides only aggregate contact
   rates by age group, the within-group heterogeneity (by activity level,
   socioeconomic status, etc.) is unknown. These bounds tell us how much this
   ignorance matters for epidemic predictions.

## Summary

| Function | What it bounds | Input needed |
|----------|---------------|--------------|
| `r0_bounds` | R₀ | Row/column sums of NGM |
| `r0_bounds_detailed_balance` | R₀ (tighter) | Row/col sums + reciprocity |
| `final_size_bounds` | Per-type τᵢ | Row or column sums + π |
| `total_final_size_bounds` | Population τ̄ | Row or column sums + π |
| `epidemic_uncertainty` | R₀ and τ̄ ranges | Sampled matrices from fiber |
