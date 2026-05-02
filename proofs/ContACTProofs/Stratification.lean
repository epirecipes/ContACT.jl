import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.BigOperators
import ContACTProofs.ContactCat

/-!
*Source: `Stratification.lean`*

# Stratification Properties

Stratification replicates a local contact structure across strata connected
by a coupling matrix:

    A_strat[(s₁,i), (s₂,j)] = C[s₁,s₂] * A[i,j]

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Stratification distributes over additive composition | ✅ |
| 2 | S commutes with ⊗ iff coupling is symmetric | ✅ |
| 3 | Symmetrisation is additive | ✅ |

-/

/-- Stratified matrix via direct product indexing. -/
noncomputable def stratifyMatrix {s a : ℕ}
    (C : Matrix (Fin s) (Fin s) ℝ)
    (A : Matrix (Fin a) (Fin a) ℝ) : Matrix (Fin s × Fin a) (Fin s × Fin a) ℝ :=
  fun p q => C p.1 q.1 * A p.2 q.2

/-- Stratification distributes over additive composition. -/
theorem stratify_distributes_over_add {s a : ℕ}
    (C : Matrix (Fin s) (Fin s) ℝ)
    (A B : Matrix (Fin a) (Fin a) ℝ) :
    stratifyMatrix C (A + B) =
    stratifyMatrix C A + stratifyMatrix C B := by
  ext p q
  simp [stratifyMatrix, Matrix.add_apply, mul_add]

/-- Symmetrisation with a population vector (local definition). -/
noncomputable def symmetriseLocal {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => (M i j * pop j + M j i * pop i) / (2 * pop j)

/-- **Key characterisation**: S(C⊗A) = C⊗S(A) at entry (s₁,i)-(s₂,j)
    iff C[s₁,s₂] = C[s₂,s₁] (given non-degenerate contact). -/
theorem symmetrise_stratify_comm_iff {s a : ℕ}
    (C : Matrix (Fin s) (Fin s) ℝ)
    (A : Matrix (Fin a) (Fin a) ℝ)
    (pop : Fin a → ℝ)
    (hpop : ∀ i, pop i > 0)
    (s₁ s₂ : Fin s) (i j : Fin a)
    (hAji : A j i * pop i ≠ 0) :
    (C s₁ s₂ * A i j * pop j + C s₂ s₁ * A j i * pop i) / (2 * pop j) =
    C s₁ s₂ * ((A i j * pop j + A j i * pop i) / (2 * pop j))
    ↔ C s₁ s₂ = C s₂ s₁ := by
  have hj_ne : pop j ≠ 0 := ne_of_gt (hpop j)
  have h2j_ne : (2 : ℝ) * pop j ≠ 0 := mul_ne_zero two_ne_zero hj_ne
  constructor
  · intro h
    -- Clear fractions: multiply both sides by 2 * pop j
    have h1 : C s₁ s₂ * A i j * pop j + C s₂ s₁ * A j i * pop i =
              C s₁ s₂ * (A i j * pop j + A j i * pop i) := by
      field_simp at h; linarith
    -- Expand RHS: C * Aij * Nj + C * Aji * Ni
    -- So: C₂₁ * Aji * Ni = C₁₂ * Aji * Ni
    have h2 : (C s₂ s₁ - C s₁ s₂) * (A j i * pop i) = 0 := by ring_nf; linarith
    rcases mul_eq_zero.mp h2 with hsub | habs
    · linarith
    · exact absurd habs hAji
  · intro hC; rw [hC]; field_simp

/-- Symmetrisation and additive composition always commute. -/
theorem symmetrise_additive_strat {n : ℕ}
    (A B : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ)
    (hpop : ∀ i, pop i > 0) :
    symmetriseLocal (A + B) pop =
    symmetriseLocal A pop + symmetriseLocal B pop := by
  ext i j
  simp [symmetriseLocal, Matrix.add_apply]
  have hj_ne : pop j ≠ 0 := ne_of_gt (hpop j)
  have h2j_ne : (2 : ℝ) * pop j ≠ 0 := mul_ne_zero two_ne_zero hj_ne
  field_simp
  ring
