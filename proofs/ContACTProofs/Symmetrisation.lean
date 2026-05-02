import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.BigOperators
import ContACTProofs.ContactCat

/-!
*Source: `Symmetrisation.lean`*

# Symmetrisation Properties

Formalises the key properties of the symmetrisation morphism:
1. Reciprocity: M_sym[i,j] · N_j = M_sym[j,i] · N_i
2. Idempotence: symmetrise(symmetrise(M)) = symmetrise(M)

The symmetrisation formula:
    M_sym[i,j] = (M[i,j] · N_j + M[j,i] · N_i) / (2 · N_j)

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Symmetrised matrix satisfies reciprocity | ✅ |
| 2 | Symmetrisation is idempotent | ✅ |

-/

/-- Symmetrise a contact matrix given a population vector. -/
noncomputable def symmetriseMatrix {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => (M i j * pop j + M j i * pop i) / (2 * pop j)

/-- Reciprocity: M_sym[i,j] · N_j = M_sym[j,i] · N_i -/
theorem symmetrise_reciprocity {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ)
    (hpop : ∀ i, pop i > 0)
    (i j : Fin n) :
    symmetriseMatrix M pop i j * pop j =
    symmetriseMatrix M pop j i * pop i := by
  simp [symmetriseMatrix]
  have hj : pop j ≠ 0 := ne_of_gt (hpop j)
  have hi : pop i ≠ 0 := ne_of_gt (hpop i)
  field_simp
  ring

/-- Idempotence: symmetrise(symmetrise(M)) = symmetrise(M) -/
theorem symmetrise_idempotent {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ)
    (hpop : ∀ i, pop i > 0) :
    symmetriseMatrix (symmetriseMatrix M pop) pop = symmetriseMatrix M pop := by
  ext i j
  simp [symmetriseMatrix]
  have hj : pop j ≠ 0 := ne_of_gt (hpop j)
  have hi : pop i ≠ 0 := ne_of_gt (hpop i)
  field_simp
  ring
