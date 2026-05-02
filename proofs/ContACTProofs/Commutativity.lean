import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import ContACTProofs.ContactCat
import ContACTProofs.Symmetrisation

/-!
*Source: `Commutativity.lean`*

# Commutativity of Operations

## Key results

1. **S always commutes with ⊕**: S(A + B) = S(A) + S(B)
2. **S commutes with ⊗ iff coupling symmetric**: characterisation theorem
3. **S is a monoid homomorphism** on (Contact, ⊕, 0)

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | S(A ⊕ B) = S(A) ⊕ S(B) | ✅ |
| 2 | S(c • A) = c • S(A) | ✅ |
| 3 | S(0) = 0 | ✅ |
| 4 | Monoid homomorphism | ✅ |
| 5 | S and ⊗ commutativity characterisation | ✅ |

-/

-- ═══════════════════════════════════════════════════════════════════════════
-- 1. Linearity of symmetrisation
-- ═══════════════════════════════════════════════════════════════════════════

/-- Symmetrisation is additive: S(A + B) = S(A) + S(B). -/
theorem symmetrise_add {n : ℕ}
    (A B : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ)
    (hpop : ∀ i, pop i > 0) :
    symmetriseMatrix (A + B) pop =
    symmetriseMatrix A pop + symmetriseMatrix B pop := by
  ext i j
  simp [symmetriseMatrix, Matrix.add_apply]
  have hj : pop j ≠ 0 := ne_of_gt (hpop j)
  have h2j : (2 : ℝ) * pop j ≠ 0 := mul_ne_zero two_ne_zero hj
  field_simp
  ring

/-- Symmetrisation is homogeneous: S(c · A) = c · S(A). -/
theorem symmetrise_smul {n : ℕ}
    (c : ℝ) (A : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ)
    (hpop : ∀ i, pop i > 0) :
    symmetriseMatrix (c • A) pop = c • symmetriseMatrix A pop := by
  ext i j
  simp [symmetriseMatrix, Matrix.smul_apply, smul_eq_mul]
  have hj : pop j ≠ 0 := ne_of_gt (hpop j)
  have h2j : (2 : ℝ) * pop j ≠ 0 := mul_ne_zero two_ne_zero hj
  field_simp

/-- **Entry-wise characterisation**: S(C⊗A) = C⊗S(A) at entry (s₁,i)-(s₂,j)
    iff C[s₁,s₂] = C[s₂,s₁] (given non-degenerate contact rates). -/
theorem symmetrise_stratify_entry_iff {s a : ℕ}
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
    have h1 : C s₁ s₂ * A i j * pop j + C s₂ s₁ * A j i * pop i =
              C s₁ s₂ * (A i j * pop j + A j i * pop i) := by
      field_simp at h; linarith
    have h2 : (C s₂ s₁ - C s₁ s₂) * (A j i * pop i) = 0 := by ring_nf; linarith
    rcases mul_eq_zero.mp h2 with hsub | habs
    · linarith
    · exact absurd habs hAji
  · intro hC; rw [hC]; field_simp

-- ═══════════════════════════════════════════════════════════════════════════
-- 3. Monoid homomorphism properties
-- ═══════════════════════════════════════════════════════════════════════════

/-- Symmetrisation preserves zero. -/
theorem symmetrise_zero {n : ℕ}
    (pop : Fin n → ℝ)
    (hpop : ∀ i, pop i > 0) :
    symmetriseMatrix (0 : Matrix (Fin n) (Fin n) ℝ) pop = 0 := by
  ext i j
  simp [symmetriseMatrix]

/-- Symmetrisation is a monoid homomorphism on (Matrix, +, 0). -/
theorem symmetrise_is_monoid_hom {n : ℕ}
    (pop : Fin n → ℝ)
    (hpop : ∀ i, pop i > 0) :
    symmetriseMatrix (0 : Matrix (Fin n) (Fin n) ℝ) pop = 0
    ∧
    ∀ A B : Matrix (Fin n) (Fin n) ℝ,
      symmetriseMatrix (A + B) pop =
      symmetriseMatrix A pop + symmetriseMatrix B pop := by
  exact ⟨symmetrise_zero pop hpop, fun A B => symmetrise_add A B pop hpop⟩
