import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import ContACTProofs.ContactCat

/-!
*Source: `ConstrainedLift.lean`*

# Constrained Lift Properties

A constrained lift reconstructs a fine matrix over a product partition P x Q
from a coarse reciprocal matrix over P and source-stratified intermediate
observations. The key categorical property is that coarsening the lift along
the projection P x Q -> P recovers the base matrix.

We formalize this in the total-contact (count) space where:
- C_fine[i,j] is the total contacts between fine groups i and j
- C_base[I,J] = sum over f(i)=I, f(j)=J of C_fine[i,j] (coarsening)
- Reciprocity: C[i,j] = C[j,i]

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Proportionate transport with equal marginals is symmetric | ✅ |
| 2 | Assembly of blocks preserves coarsening | ✅ |
| 3 | Structural zeros (fiber) are preserved by assembly | ✅ |

-/

/-- Proportionate transport: C[i,j] = r[i] * c[j] / T.
    The unique matrix with given row and column marginals under
    proportionate mixing. -/
noncomputable def proportionateTransport {n m : ℕ}
    (r : Fin n → ℝ) (c : Fin m → ℝ) (T : ℝ) : Matrix (Fin n) (Fin m) ℝ :=
  fun i j => if T = 0 then 0 else r i * c j / T

/-- Proportionate transport with equal row and column marginals is symmetric. -/
theorem proportionate_symmetric {n : ℕ}
    (r : Fin n → ℝ)
    (T : ℝ)
    (hT : T ≠ 0) :
    Matrix.transpose (proportionateTransport r r T) =
    proportionateTransport r r T := by
  ext i j
  simp only [Matrix.transpose_apply, proportionateTransport, if_neg hT]
  ring

/-- A partition map f : Fin m -> Fin n induces coarsening by summing over fibers. -/
noncomputable def coarsenCounts {m n : ℕ}
    (C : Matrix (Fin m) (Fin m) ℝ)
    (f : Fin m → Fin n) : Matrix (Fin n) (Fin n) ℝ :=
  fun I J => ∑ i ∈ Finset.univ.filter (fun i => f i = I),
             ∑ j ∈ Finset.univ.filter (fun j => f j = J),
             C i j

/-- If we construct C_fine by placing blocks for each base pair (I,J) such that
    the block sum equals blockTotal I J, then coarsening recovers blockTotal. -/
theorem coarsen_block_assembly {m n : ℕ}
    (f : Fin m → Fin n)
    (blockTotal : Fin n → Fin n → ℝ)
    (C : Matrix (Fin m) (Fin m) ℝ)
    (hblock : ∀ I J,
      ∑ i ∈ Finset.univ.filter (fun i => f i = I),
      ∑ j ∈ Finset.univ.filter (fun j => f j = J),
      C i j = blockTotal I J) :
    coarsenCounts C f = fun I J => blockTotal I J := by
  ext I J
  simp only [coarsenCounts]
  exact hblock I J

/-- Structural zeros (fiber version): if C[i, .] = 0 and C[., j] = 0 for ALL
    elements in the fiber of a base group I0, then the coarsened matrix has
    zero row and column at I0. -/
theorem structural_zeros_fiber {m n : ℕ}
    (f : Fin m → Fin n)
    (C : Matrix (Fin m) (Fin m) ℝ)
    (I0 : Fin n)
    (hrow : ∀ i, f i = I0 → ∀ j, C i j = 0)
    (hcol : ∀ j, f j = I0 → ∀ i, C i j = 0) :
    (∀ J, coarsenCounts C f I0 J = 0) ∧
    (∀ I, coarsenCounts C f I I0 = 0) := by
  constructor
  · intro J
    simp only [coarsenCounts]
    apply Finset.sum_eq_zero
    intro i hi
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hi
    apply Finset.sum_eq_zero
    intro j _
    exact hrow i hi j
  · intro I
    simp only [coarsenCounts]
    apply Finset.sum_eq_zero
    intro i _
    apply Finset.sum_eq_zero
    intro j hj
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hj
    exact hcol j hj i

/-- If C is symmetric and we coarsen it, the result is symmetric. -/
theorem coarsen_preserves_symmetry {m n : ℕ}
    (C : Matrix (Fin m) (Fin m) ℝ)
    (f : Fin m → Fin n)
    (hsym : ∀ i j, C i j = C j i) :
    ∀ I J, coarsenCounts C f I J = coarsenCounts C f J I := by
  intro I J
  simp only [coarsenCounts]
  conv_rhs => rw [Finset.sum_comm]
  congr 1
  ext i
  congr 1
  ext j
  exact hsym i j
