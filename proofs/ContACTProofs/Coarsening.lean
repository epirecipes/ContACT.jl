import Mathlib.CategoryTheory.Category.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import ContACTProofs.ContactCat

/-!
*Source: `Coarsening.lean`*

# Coarsening Functoriality

Coarsening is functorial: given age-group maps f and g, we have
coarsen(g ∘ f) = coarsen(g) ∘ coarsen(f).

We use the "total contacts" formulation where coarsening sums matrix entries
over fibers — this avoids division and makes functoriality follow from
basic sum manipulations.

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Coarsening along identity = identity | ✅ |
| 2 | Coarsening preserves total contacts | ✅ |
| 3 | Coarsening is functorial | ✅ |

-/

/-- Coarsening a matrix by summing entries over fibers of f. -/
noncomputable def coarsenTotal {m n : ℕ}
    (C : Matrix (Fin m) (Fin m) ℝ)
    (f : Fin m → Fin n) : Matrix (Fin n) (Fin n) ℝ :=
  fun I J => ∑ i ∈ Finset.univ.filter (fun i => f i = I),
             ∑ j ∈ Finset.univ.filter (fun j => f j = J),
             C i j

/-- Coarsening along the identity is the identity. -/
theorem coarsen_id {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ) :
    coarsenTotal C id = C := by
  ext I J
  simp only [coarsenTotal, Function.id_def]
  rw [show Finset.univ.filter (fun i => i = I) = {I} from by ext x; simp]
  rw [show Finset.univ.filter (fun j => j = J) = {J} from by ext x; simp]
  simp

/-- Total contacts are preserved under coarsening:
    ∑_{I,J} coarsen(C)[I,J] = ∑_{i,j} C[i,j].
    The fibers of f partition Fin m, so summing over coarse groups and
    fibers is equivalent to summing over all fine indices. -/
theorem coarsen_preserves_total {m n : ℕ}
    (C : Matrix (Fin m) (Fin m) ℝ)
    (f : Fin m → Fin n) :
    ∑ I : Fin n, ∑ J : Fin n, coarsenTotal C f I J =
    ∑ i : Fin m, ∑ j : Fin m, C i j := by
  simp only [coarsenTotal]
  -- LHS: ∑ I, ∑ J, ∑ i ∈ fib(f, I), ∑ j ∈ fib(f, J), C i j
  -- Swap inner two: ∑ I, ∑ i ∈ fib(f,I), ∑ J, ∑ j ∈ fib(f, J), C i j
  conv_lhs =>
    arg 2; ext I
    rw [Finset.sum_comm]
  -- Apply sum_fiberwise to inner (J, j): ∑ J, ∑ j ∈ fib(f, J), C i j = ∑ j, C i j
  conv_lhs =>
    arg 2; ext I; arg 2; ext i
    rw [Finset.sum_fiberwise Finset.univ f (fun j => C i j)]
  -- Now: ∑ I, ∑ i ∈ fib(f, I), ∑ j, C i j = ∑ i, ∑ j, C i j
  exact Finset.sum_fiberwise Finset.univ f (fun i => ∑ j, C i j)

/-- Helper: sum over fibers of (g ∘ f) decomposes as nested sum over fibers of f then g. -/
private lemma sum_fiber_comp {l m n : ℕ}
    (f : Fin l → Fin m) (g : Fin m → Fin n) (K : Fin n) (h : Fin l → ℝ) :
    ∑ i ∈ Finset.univ.filter (fun i => g (f i) = K), h i =
    ∑ I ∈ Finset.univ.filter (fun I => g I = K),
      ∑ i ∈ Finset.univ.filter (fun i => f i = I), h i := by
  rw [show ∑ I ∈ Finset.univ.filter (fun I => g I = K),
        ∑ i ∈ Finset.univ.filter (fun i => f i = I), h i =
      ∑ i ∈ Finset.univ.filter (fun i => g (f i) = K), h i from ?_]
  -- Use sum_fiberwise on the right side after relating the filters
  rw [show Finset.univ.filter (fun i : Fin l => g (f i) = K) =
        (Finset.univ.filter (fun I : Fin m => g I = K)).biUnion
          (fun I => Finset.univ.filter (fun i : Fin l => f i = I)) from ?_]
  · rw [Finset.sum_biUnion]
    intros I _ J _ hIJ
    simp only [Finset.disjoint_filter]
    intros i _ hfi
    intro hfi'
    exact hIJ (hfi.symm.trans hfi')
  · ext i
    simp only [Finset.mem_filter, Finset.mem_univ, Finset.mem_biUnion, true_and]
    constructor
    · intro hi; exact ⟨f i, hi, rfl⟩
    · rintro ⟨I, hI, hfi⟩; rw [hfi]; exact hI

/-- **Functoriality**: coarsening along a composition equals composing coarsenings.
    coarsenTotal C (g ∘ f) = coarsenTotal (coarsenTotal C f) g -/
theorem coarsen_functorial {l m n : ℕ}
    (C : Matrix (Fin l) (Fin l) ℝ)
    (f : Fin l → Fin m) (g : Fin m → Fin n) :
    coarsenTotal C (g ∘ f) =
    coarsenTotal (coarsenTotal C f) g := by
  ext K L
  simp only [coarsenTotal, Function.comp]
  -- LHS: ∑ i ∈ fib(g∘f, K), ∑ j ∈ fib(g∘f, L), C i j
  -- RHS: ∑ I ∈ fib(g, K), ∑ J ∈ fib(g, L),
  --        ∑ i ∈ fib(f, I), ∑ j ∈ fib(f, J), C i j
  -- Apply sum_fiber_comp on the i sum (LHS)
  rw [sum_fiber_comp f g K (fun i => ∑ j ∈ Finset.univ.filter (fun j => g (f j) = L), C i j)]
  -- LHS now: ∑ I ∈ fib(g, K), ∑ i ∈ fib(f, I), ∑ j ∈ fib(g∘f, L), C i j
  apply Finset.sum_congr rfl; intro I _
  -- Goal: ∑ i ∈ fib(f, I), ∑ j ∈ fib(g∘f, L), C i j
  --     = ∑ J ∈ fib(g, L), ∑ i ∈ fib(f, I), ∑ j ∈ fib(f, J), C i j
  -- Apply sum_fiber_comp to the inner j sum
  conv_lhs =>
    arg 2; ext i
    rw [sum_fiber_comp f g L (fun j => C i j)]
  -- LHS now: ∑ i ∈ fib(f, I), ∑ J ∈ fib(g, L), ∑ j ∈ fib(f, J), C i j
  -- Swap i and J
  exact Finset.sum_comm
