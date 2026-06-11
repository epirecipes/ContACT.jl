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

-- ═══════════════════════════════════════════════════════════════════════════
-- Population-weighted coarsening (the operation actually implemented in
-- `coarsen` in `src/coarsening.jl`)
-- ═══════════════════════════════════════════════════════════════════════════

/-!
The Julia `coarsen` operates on a *mean-contacts* matrix `M` with a population
vector `N`:

    M_c[I,J] = (∑_{j∈f⁻¹J} N_j · ∑_{i∈f⁻¹I} M[i,j]) / (∑_{j∈f⁻¹J} N_j).

This is the population-weighted average — not the plain fiber sum `coarsenTotal`.
We model it here and prove that it is functorial and preserves total contacts,
so the README claims hold for the operation the code ships, not only for the
simplified count model. The bridge is that coarsening total contacts
`Counts[i,j] = M[i,j]·N_j` over fibers is exactly `coarsenTotal Counts`, and the
mean-coarsening is that divided by the coarse population.
-/

/-- Population of a coarse group: the sum of fine populations in its fiber. -/
noncomputable def coarsenPop {m n : ℕ} (N : Fin m → ℝ) (f : Fin m → Fin n) :
    Fin n → ℝ :=
  fun J => ∑ j ∈ Finset.univ.filter (fun j => f j = J), N j

/-- Population-weighted coarsening of a mean-contacts matrix, matching
    `coarsen` in `src/coarsening.jl`. -/
noncomputable def weightedCoarsen {m n : ℕ}
    (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun I J => coarsenTotal (fun i j => M i j * N j) f I J / coarsenPop N f J

theorem weightedCoarsen_apply {m n : ℕ}
    (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n) (I J : Fin n) :
    weightedCoarsen M N f I J =
    coarsenTotal (fun i j => M i j * N j) f I J / coarsenPop N f J := rfl

/-- Coarse total contacts equal the summed fine total contacts over the fiber:
    `M_c[I,J] · N_J = ∑_{i∈f⁻¹I, j∈f⁻¹J} M[i,j]·N_j`.
    Non-negativity of populations means an empty coarse group (`N_J = 0`) carries
    no contacts, so the identity also holds there. -/
theorem weightedCoarsen_mul_pop {m n : ℕ}
    (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) (I J : Fin n) :
    weightedCoarsen M N f I J * coarsenPop N f J =
    coarsenTotal (fun i j => M i j * N j) f I J := by
  unfold weightedCoarsen
  rcases eq_or_ne (coarsenPop N f J) 0 with hz | hz
  · rw [hz, mul_zero]
    symm
    simp only [coarsenTotal]
    have hzero : ∀ j ∈ Finset.univ.filter (fun j => f j = J), N j = 0 :=
      (Finset.sum_eq_zero_iff_of_nonneg (fun j _ => hN j)).mp (by simpa [coarsenPop] using hz)
    apply Finset.sum_eq_zero
    intro i _
    apply Finset.sum_eq_zero
    intro j hj
    rw [hzero j hj, mul_zero]
  · field_simp

/-- **Total contacts are preserved** by the population-weighted coarsening:
    `∑_{I,J} M_c[I,J]·N_J = ∑_{i,j} M[i,j]·N_j`. This is the README's
    "coarsening preserves total contacts" for the shipped `coarsen`. -/
theorem weightedCoarsen_preserves_total {m n : ℕ}
    (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) :
    (∑ I : Fin n, ∑ J : Fin n, weightedCoarsen M N f I J * coarsenPop N f J) =
    ∑ i : Fin m, ∑ j : Fin m, M i j * N j := by
  have h : ∀ I J : Fin n,
      weightedCoarsen M N f I J * coarsenPop N f J =
      coarsenTotal (fun i j => M i j * N j) f I J :=
    fun I J => weightedCoarsen_mul_pop M N f hN I J
  simp_rw [h]
  exact coarsen_preserves_total (fun i j => M i j * N j) f

/-- Coarse populations telescope along a composition: aggregating fine
    populations to the intermediate then the coarse level equals aggregating
    directly. -/
theorem coarsenPop_comp {l m n : ℕ}
    (N : Fin l → ℝ) (f : Fin l → Fin m) (g : Fin m → Fin n) (K : Fin n) :
    coarsenPop (coarsenPop N f) g K = coarsenPop N (g ∘ f) K := by
  simp only [coarsenPop, Function.comp_apply]
  exact (sum_fiber_comp f g K N).symm

/-- **Functoriality of the population-weighted coarsening**:
    `weightedCoarsen M N (g ∘ f) = weightedCoarsen (weightedCoarsen M N f) (coarsenPop N f) g`.
    The intermediate population is the aggregated `coarsenPop N f`. This is the
    README's coarsening functoriality for the operation the code ships. -/
theorem weightedCoarsen_functorial {l m n : ℕ}
    (M : Matrix (Fin l) (Fin l) ℝ) (N : Fin l → ℝ)
    (f : Fin l → Fin m) (g : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) :
    weightedCoarsen M N (g ∘ f) =
    weightedCoarsen (weightedCoarsen M N f) (coarsenPop N f) g := by
  ext K L
  rw [weightedCoarsen_apply, weightedCoarsen_apply, coarsenPop_comp]
  have hpt : (fun (I J : Fin m) => weightedCoarsen M N f I J * coarsenPop N f J)
      = coarsenTotal (fun i j => M i j * N j) f := by
    funext I J
    exact weightedCoarsen_mul_pop M N f hN I J
  rw [hpt, ← coarsen_functorial]

-- ═══════════════════════════════════════════════════════════════════════════
-- The coarsening action of category morphisms (links `ContactCat` to coarsening)
-- ═══════════════════════════════════════════════════════════════════════════

/-- A contact-matrix morphism acts on the (count-space) contact matrix of its
    domain by coarsening along the underlying age-group map. This gives the
    morphisms of `ContactCat` a genuine action on contact-matrix data, beyond the
    bare `Fin`-map. -/
noncomputable def ContactMatHom.coarsenAction {A B : ContactMatObj}
    (φ : ContactMatHom A B) (C : Matrix (Fin A.n) (Fin A.n) ℝ) :
    Matrix (Fin B.n) (Fin B.n) ℝ :=
  coarsenTotal C φ.map

/-- The identity morphism acts as the identity on contact data. -/
theorem ContactMatHom.coarsenAction_id {A : ContactMatObj}
    (C : Matrix (Fin A.n) (Fin A.n) ℝ) :
    (ContactMatHom.id A).coarsenAction C = C :=
  coarsen_id C

/-- The action respects composition: coarsening along `g ∘ f` equals coarsening
    along `f` then `g`. Together with `coarsenAction_id` this makes the action a
    functor on **Contact**, so "contact matrices form a category" carries
    contact-structure content rather than just `Fin`-map composition. -/
theorem ContactMatHom.coarsenAction_comp {A B C : ContactMatObj}
    (g : ContactMatHom B C) (f : ContactMatHom A B)
    (M : Matrix (Fin A.n) (Fin A.n) ℝ) :
    (g.comp f).coarsenAction M = g.coarsenAction (f.coarsenAction M) := by
  simp only [ContactMatHom.coarsenAction, ContactMatHom.comp]
  exact coarsen_functorial M f.map g.map
