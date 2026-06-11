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

-- ═══════════════════════════════════════════════════════════════════════════
-- Zero-population branch (matching `symmetrise` in `src/symmetrise.jl`)
-- ═══════════════════════════════════════════════════════════════════════════

/-!
The theorems above assume every group has strictly positive population. The Julia
`symmetrise` instead handles empty groups explicitly: when `N_j = 0` it sets
`M_sym[i,j] = 0`, and it *requires* the corresponding total contacts to vanish
(throwing otherwise). We model that branch with `symmetriseSafe` and prove
reciprocity and idempotence under the same consistency condition the code
enforces — namely that an empty group carries no contacts:
`pop a = 0 → M a b · pop b = 0`. This covers the zero-population path the
strictly-positive theorems exclude.
-/

/-- Symmetrisation matching the Julia code's zero-population branch: empty groups
    (`pop j = 0`) get a zero column. -/
noncomputable def symmetriseSafe {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => if pop j = 0 then 0 else (M i j * pop j + M j i * pop i) / (2 * pop j)

theorem symmetriseSafe_apply {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop : Fin n → ℝ) (i j : Fin n) :
    symmetriseSafe M pop i j =
    if pop j = 0 then 0 else (M i j * pop j + M j i * pop i) / (2 * pop j) := rfl

/-- On nonempty groups `symmetriseSafe` agrees with the plain `symmetriseMatrix`. -/
theorem symmetriseSafe_eq_of_pos {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop : Fin n → ℝ) (i j : Fin n) (h : pop j ≠ 0) :
    symmetriseSafe M pop i j = symmetriseMatrix M pop i j := by
  rw [symmetriseSafe_apply, if_neg h]; rfl

/-- Total contacts of the symmetrised matrix are the symmetric average
    `(M[i,j]·N_j + M[j,i]·N_i)/2`, in both the populated and empty-group cases.
    In the empty case the consistency condition forces the total to vanish. -/
theorem symmetriseSafe_mul_pop {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop : Fin n → ℝ)
    (hcons : ∀ a b, pop a = 0 → M a b * pop b = 0)
    (i j : Fin n) :
    symmetriseSafe M pop i j * pop j = (M i j * pop j + M j i * pop i) / 2 := by
  rw [symmetriseSafe_apply]
  by_cases hj : pop j = 0
  · rw [if_pos hj, hj, hcons j i hj]; ring
  · rw [if_neg hj]
    field_simp

/-- **Reciprocity** including empty groups:
    `M_sym[i,j]·N_j = M_sym[j,i]·N_i`. -/
theorem symmetriseSafe_reciprocity {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop : Fin n → ℝ)
    (hcons : ∀ a b, pop a = 0 → M a b * pop b = 0)
    (i j : Fin n) :
    symmetriseSafe M pop i j * pop j = symmetriseSafe M pop j i * pop i := by
  rw [symmetriseSafe_mul_pop M pop hcons i j, symmetriseSafe_mul_pop M pop hcons j i]
  ring

/-- **Idempotence** including empty groups, matching the shipped `symmetrise`:
    `symmetrise(symmetrise(M)) = symmetrise(M)` for all non-negative populations
    consistent with the code's zero-total requirement. -/
theorem symmetriseSafe_idempotent {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop : Fin n → ℝ)
    (hcons : ∀ a b, pop a = 0 → M a b * pop b = 0) :
    symmetriseSafe (symmetriseSafe M pop) pop = symmetriseSafe M pop := by
  ext i j
  rw [symmetriseSafe_apply (symmetriseSafe M pop) pop i j]
  by_cases hj : pop j = 0
  · rw [if_pos hj, symmetriseSafe_apply M pop i j, if_pos hj]
  · rw [if_neg hj, symmetriseSafe_mul_pop M pop hcons i j,
        symmetriseSafe_mul_pop M pop hcons j i, symmetriseSafe_apply M pop i j, if_neg hj]
    ring
