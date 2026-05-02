import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import ContACTProofs.ContactCat

/-!
*Source: `Composition.lean`*

# Additive Composition Laws

Contact matrices from different settings (home, work, school, other) compose
additively. This file proves the monoid laws:
1. Associativity: (A ⊕ B) ⊕ C = A ⊕ (B ⊕ C)
2. Commutativity: A ⊕ B = B ⊕ A
3. Identity: A ⊕ 0 = A

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Additive composition is associative | ✅ |
| 2 | Additive composition is commutative | ✅ |
| 3 | Zero matrix is the identity element | ✅ |

-/

/-- Additive composition of contact matrices (setting combination). -/
def composeContact {n : ℕ}
    (A B : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  A + B

/-- Associativity of additive composition. -/
theorem compose_assoc {n : ℕ}
    (A B C : Matrix (Fin n) (Fin n) ℝ) :
    composeContact (composeContact A B) C =
    composeContact A (composeContact B C) := by
  simp [composeContact, add_assoc]

/-- Commutativity of additive composition. -/
theorem compose_comm {n : ℕ}
    (A B : Matrix (Fin n) (Fin n) ℝ) :
    composeContact A B = composeContact B A := by
  simp [composeContact, add_comm]

/-- Zero matrix is the identity for composition. -/
theorem compose_zero_left {n : ℕ}
    (A : Matrix (Fin n) (Fin n) ℝ) :
    composeContact 0 A = A := by
  simp [composeContact]

theorem compose_zero_right {n : ℕ}
    (A : Matrix (Fin n) (Fin n) ℝ) :
    composeContact A 0 = A := by
  simp [composeContact]
