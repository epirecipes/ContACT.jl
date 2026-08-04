import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Algebra.BigOperators.Field
import ContACTProofs.ContactCat
import ContACTProofs.Symmetrisation

/-!
*Source: `Reprojection.lean`*

# Reprojection Properties

`reproject` (`src/reprojection.jl`) is `symmetriseMatrix` evaluated at a target
population `target` in place of a matrix's own population. This file proves:
1. Reciprocity: M'[i,j]·N'_j = M'[j,i]·N'_i, including empty target groups
2. Total contacts preserved: Σ M'[i,j]·N'_j = Σ M[i,j]·N'_j, including empty
   target groups
3. Identity on matrices already reciprocal at `target`
4. Not functorial in the population

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Reprojected matrix satisfies reciprocity at the target population | ✅ |
| 2 | Reciprocity including empty target groups | ✅ |
| 3 | Total contacts at the target population are preserved | ✅ |
| 4 | Total contacts preserved including empty target groups | ✅ |
| 5 | Identity on matrices already reciprocal at the target population | ✅ |
| 6 | Reprojection is not functorial in the population, even from a matrix already reciprocal at its own population | ✅ |

-/

/-- Reciprocity: M'[i,j]·N'_j = M'[j,i]·N'_i. -/
theorem reproject_reciprocity {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (target : Fin n → ℝ)
    (htarget : ∀ i, target i > 0) (i j : Fin n) :
    symmetriseMatrix M target i j * target j =
    symmetriseMatrix M target j i * target i :=
  symmetrise_reciprocity M target htarget i j

/-- Reciprocity including empty target groups, matching `_symmetrise_at`
(`src/symmetrise.jl`, shared by `symmetrise` and `reproject`). -/
theorem reproject_reciprocity_safe {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (target : Fin n → ℝ)
    (hcons : ∀ a b, target a = 0 → M a b * target b = 0) (i j : Fin n) :
    symmetriseSafe M target i j * target j =
    symmetriseSafe M target j i * target i :=
  symmetriseSafe_reciprocity M target hcons i j

/-- Total contacts at the target population are preserved:
Σ_{i,j} M'[i,j]·N'_j = Σ_{i,j} M[i,j]·N'_j. -/
theorem reproject_total_contacts {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (target : Fin n → ℝ) (htarget : ∀ j, target j ≠ 0) :
    (∑ i : Fin n, ∑ j : Fin n, symmetriseMatrix M target i j * target j) =
    ∑ i : Fin n, ∑ j : Fin n, M i j * target j := by
  have hentry : ∀ i j : Fin n,
      symmetriseMatrix M target i j * target j = (M i j * target j + M j i * target i) / 2 := by
    intro i j
    have hj := htarget j
    simp only [symmetriseMatrix]
    field_simp
  simp_rw [hentry, add_div, Finset.sum_add_distrib]
  rw [Finset.sum_comm (f := fun i j => M j i * target i / 2)]
  simp_rw [← Finset.sum_div]
  ring

/-- Total contacts at the target population are preserved, including empty target
groups (`_symmetrise_at`'s zero-population branch). -/
theorem reproject_total_contacts_safe {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (target : Fin n → ℝ)
    (hcons : ∀ a b, target a = 0 → M a b * target b = 0) :
    (∑ i : Fin n, ∑ j : Fin n, symmetriseSafe M target i j * target j) =
    ∑ i : Fin n, ∑ j : Fin n, M i j * target j := by
  simp_rw [symmetriseSafe_mul_pop M target hcons, add_div, Finset.sum_add_distrib]
  rw [Finset.sum_comm (f := fun i j => M j i * target i / 2)]
  simp_rw [← Finset.sum_div]
  ring

/-- Identity on matrices already reciprocal at `target`. -/
theorem reproject_id_of_reciprocal {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (target : Fin n → ℝ) (htarget : ∀ i, target i ≠ 0)
    (hrec : ∀ i j, M i j * target j = M j i * target i) :
    symmetriseMatrix M target = M := by
  ext i j
  have hj := htarget j
  simp only [symmetriseMatrix]
  rw [← hrec i j]
  field_simp
  ring

/-- Not functorial in the population, even starting from a matrix already
reciprocal at its own population. Witnesses: `M = !![0,1;1,0]`, reciprocal at
`N = (1,1)`; reprojecting through `N' = (1,2)` then onto `N'' = (1,1)`
disagrees with reprojecting directly onto `N''`. Two-step gives entry
`(0,1) = 9/8`; direct gives `1`. -/
theorem reproject_not_functorial :
    ∃ (M : Matrix (Fin 2) (Fin 2) ℝ) (N N' N'' : Fin 2 → ℝ),
      (∀ i, 0 < N i) ∧ (∀ i, 0 < N' i) ∧ (∀ i, 0 < N'' i) ∧
      (∀ i j, M i j * N j = M j i * N i) ∧
      symmetriseMatrix (symmetriseMatrix M N') N'' ≠ symmetriseMatrix M N'' := by
  refine ⟨!![0, 1; 1, 0], ![1, 1], ![1, 2], ![1, 1], ?_, ?_, ?_, ?_, ?_⟩
  · intro i; fin_cases i <;> norm_num
  · intro i; fin_cases i <;> norm_num
  · intro i; fin_cases i <;> norm_num
  · intro i j; fin_cases i <;> fin_cases j <;> norm_num
  · intro h
    have h01 := congrFun (congrFun h 0) 1
    simp [symmetriseMatrix] at h01
    norm_num at h01
