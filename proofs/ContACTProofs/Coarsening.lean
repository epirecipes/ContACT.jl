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
| 4 | Population-weighted coarsening (`weightedCoarsen`, matching the shipped `coarsen`) preserves total contacts | ✅ |
| 5 | Population-weighted coarsening is functorial | ✅ |
| 6 | `ContactCat` morphisms act on contact data via coarsening, functorially | ✅ |
| 7 | Activity-vector build/coarsen equivariance: the reciprocal proportionate-mixing lift `a_i·N_i·a_j/D` commutes with coarsening for any fibre map (`weightedCoarsen_proportionate`, and `_nonneg` with no coarse-population condition) | ✅ |

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
-- ═══════════════════════════════════════════════════════════════════════════
-- Build/coarsen equivariance for activity-vector lifts
-- ═══════════════════════════════════════════════════════════════════════════

/-!
An activity vector `a` (mean contacts per person, per group) can be turned
into a ContactMatrix, and a ContactMatrix can be coarsened along a fibre map
that merges fine groups into coarse ones. These two operations should commute:
coarsening the activity vector and then building should give the same matrix
as building and then coarsening. If they don't, a model fitted at one spatial
or demographic resolution can't be interpreted at another.

Three findings pinned here:
  1. The reciprocal lift `M[i,j] = a_i·N_i·a_j / D` (`D = Σ_k a_k·N_k`) commutes
     for any fibre map, contiguous or not, needing only that each coarse group
     has nonzero total population. Both paths reduce to `S_I·S_J/(D·N_J)` with
     `S_I = Σ_{i∈I} a_i·N_i`. This is `weightedCoarsen_proportionate` below.
  2. The naive lift `M[i,j] = a_i·a_j/D` (missing the `N_i` factor) is not a
     valid MeanContacts ContactMatrix: MeanContacts reciprocity requires
     symmetry of *total* contacts, `M[i,j]·N_j == M[j,i]·N_i`, which fails
     whenever group populations differ. So it is a trap under MeanContacts,
     not a counterexample to (1). The factor is required by ContACT's
     convention that column = participant (per-capita) and row = contactee
     (extensive).
  3. The Britton-Ball assortative/disassortative kernels allocate contacts by
     walking activity strata in order, so they commute with coarsening only
     when the fibre map is order-preserving (each coarse group is a contiguous
     run of fine strata). Merging non-adjacent strata changes the order the
     algorithm depends on and the mismatch is O(1) in mean-contacts units, not
     floating-point noise. Normal usage merges adjacent strata, so this is a
     documented gap rather than a blocker — see the Julia regression test in
     `test/runtests.jl` for (2) and (3); only (1) is formalised here.
-/

/-- Total (population-weighted) activity `D = Σ_k a_k·N_k`, the normaliser of
    proportionate mixing. -/
noncomputable def totalActivity {n : ℕ} (a N : Fin n → ℝ) : ℝ :=
  ∑ k, a k * N k

/-- The reciprocal proportionate-mixing lift of a per-capita activity vector:
    `M[i,j] = a_i·N_i·a_j / D`. -/
noncomputable def proportionateMatrix {n : ℕ} (a N : Fin n → ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun i j => a i * N i * a j / totalActivity a N

/-- Population-weighted activity mass in a fibre: `S_I = Σ_{i∈f⁻¹I} a_i·N_i`. -/
noncomputable def coarsenActivityMass {m n : ℕ} (a N : Fin m → ℝ) (f : Fin m → Fin n) :
    Fin n → ℝ :=
  fun I => ∑ i ∈ Finset.univ.filter (fun i => f i = I), a i * N i

/-- Population-weighted mean pushforward of a per-capita activity vector along a
    fibre map: `a_I = S_I / N_I`. This is the natural pushforward of an
    *intensive* quantity, as opposed to `coarsenPop`, which sums the *extensive*
    population. -/
noncomputable def coarsenActivity {m n : ℕ} (a N : Fin m → ℝ) (f : Fin m → Fin n) :
    Fin n → ℝ :=
  fun I => coarsenActivityMass a N f I / coarsenPop N f I

theorem coarsenActivity_apply {m n : ℕ} (a N : Fin m → ℝ) (f : Fin m → Fin n) (I : Fin n) :
    coarsenActivity a N f I = coarsenActivityMass a N f I / coarsenPop N f I := rfl

/-- `coarsenTotal` distributes over division by a constant. -/
theorem coarsenTotal_div_const {m n : ℕ} (C : Matrix (Fin m) (Fin m) ℝ) (D : ℝ)
    (f : Fin m → Fin n) (I J : Fin n) :
    coarsenTotal (fun i j => C i j / D) f I J = coarsenTotal C f I J / D := by
  simp only [coarsenTotal, div_eq_mul_inv, ← Finset.sum_mul]

/-- `coarsenTotal` of an outer product factors into the product of the two
    fibre sums. -/
theorem coarsenTotal_outer {m n : ℕ} (u v : Fin m → ℝ) (f : Fin m → Fin n) (I J : Fin n) :
    coarsenTotal (fun i j => u i * v j) f I J =
    (∑ i ∈ Finset.univ.filter (fun i => f i = I), u i) *
    (∑ j ∈ Finset.univ.filter (fun j => f j = J), v j) := by
  show (∑ i ∈ Finset.univ.filter (fun i => f i = I),
          ∑ j ∈ Finset.univ.filter (fun j => f j = J), u i * v j) = _
  rw [Finset.sum_mul]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Finset.mul_sum]

/-- Combined form used by the main theorem. -/
theorem coarsenTotal_outer_div {m n : ℕ} (u v : Fin m → ℝ) (D : ℝ) (f : Fin m → Fin n) (I J : Fin n) :
    coarsenTotal (fun i j => u i * v j / D) f I J =
    (∑ i ∈ Finset.univ.filter (fun i => f i = I), u i) *
    (∑ j ∈ Finset.univ.filter (fun j => f j = J), v j) / D := by
  rw [coarsenTotal_div_const (fun i j => u i * v j) D f I J, coarsenTotal_outer]

/-- Total activity is preserved by the population-weighted pushforward:
    `D` computed from the coarse `(coarsenActivity a N f, coarsenPop N f)`
    equals `D` computed from the fine `(a, N)`. -/
theorem totalActivity_coarsen {m n : ℕ} (a N : Fin m → ℝ) (f : Fin m → Fin n)
    (hNJ : ∀ J : Fin n, coarsenPop N f J ≠ 0) :
    totalActivity (coarsenActivity a N f) (coarsenPop N f) = totalActivity a N := by
  have hcancel : ∀ I : Fin n,
      coarsenActivity a N f I * coarsenPop N f I = coarsenActivityMass a N f I := by
    intro I
    rw [coarsenActivity_apply, div_mul_cancel₀ _ (hNJ I)]
  show ∑ I : Fin n, coarsenActivity a N f I * coarsenPop N f I = ∑ i : Fin m, a i * N i
  simp_rw [hcancel]
  exact Finset.sum_fiberwise Finset.univ f (fun i => a i * N i)

/-- The cancellation `a_I · N_I = S_I` underlying both coarsening-equivariance
    results, without assuming the coarse population is nonzero. Populations are
    nonnegative in every `ContactMatrix`, so a fibre summing to zero forces each of
    its members to zero; the activity mass then vanishes too and both sides are `0`.
    This is what lets `hNJ` be discharged rather than merely documented. -/
theorem coarsenActivity_mul_pop {m n : ℕ} (a N : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) (I : Fin n) :
    coarsenActivity a N f I * coarsenPop N f I = coarsenActivityMass a N f I := by
  rcases eq_or_ne (coarsenPop N f I) 0 with hz | hz
  · rw [hz, mul_zero]
    symm
    simp only [coarsenPop] at hz
    have hzero : ∀ i ∈ Finset.univ.filter (fun i => f i = I), N i = 0 :=
      (Finset.sum_eq_zero_iff_of_nonneg (fun i _ => hN i)).mp hz
    simp only [coarsenActivityMass]
    exact Finset.sum_eq_zero (fun i hi => by rw [hzero i hi, mul_zero])
  · rw [coarsenActivity_apply, div_mul_cancel₀ _ hz]

/-- Total activity is preserved by coarsening, on nonnegative populations and with
    no condition on the coarse populations. -/
theorem totalActivity_coarsen_nonneg {m n : ℕ} (a N : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) :
    totalActivity (coarsenActivity a N f) (coarsenPop N f) = totalActivity a N := by
  show ∑ I : Fin n, coarsenActivity a N f I * coarsenPop N f I = ∑ i : Fin m, a i * N i
  simp_rw [coarsenActivity_mul_pop a N f hN]
  exact Finset.sum_fiberwise Finset.univ f (fun i => a i * N i)

/-- **Build/coarsen equivariance on the domain the package actually uses.**
    Identical conclusion to `weightedCoarsen_proportionate`, with the
    coarse-population hypothesis `hNJ` replaced by nonnegativity of the fine
    populations — which every `ContactMatrix` satisfies by construction
    (`src/types.jl` rejects negative entries). On that domain the equivariance holds
    with no side condition beyond the populations being populations, for any fibre
    map, surjective or not. -/
theorem weightedCoarsen_proportionate_nonneg {m n : ℕ}
    (a N : Fin m → ℝ) (f : Fin m → Fin n) (hN : ∀ j, 0 ≤ N j) :
    weightedCoarsen (proportionateMatrix a N) N f =
    proportionateMatrix (coarsenActivity a N f) (coarsenPop N f) := by
  have hD := totalActivity_coarsen_nonneg a N f hN
  ext I J
  have hfun : (fun i j => proportionateMatrix a N i j * N j)
      = (fun i j => (a i * N i) * (a j * N j) / totalActivity a N) := by
    funext i j
    unfold proportionateMatrix
    ring
  show coarsenTotal (fun i j => proportionateMatrix a N i j * N j) f I J / coarsenPop N f J
     = proportionateMatrix (coarsenActivity a N f) (coarsenPop N f) I J
  rw [hfun, coarsenTotal_outer_div]
  show coarsenActivityMass a N f I * coarsenActivityMass a N f J / totalActivity a N
       / coarsenPop N f J
     = coarsenActivity a N f I * coarsenPop N f I * coarsenActivity a N f J
       / totalActivity (coarsenActivity a N f) (coarsenPop N f)
  rw [hD, coarsenActivity_mul_pop a N f hN I, coarsenActivity_apply]
  ring

/-- Coarsening commutes with the reciprocal proportionate-mixing lift:
    coarsening an activity vector and then building a matrix agrees with
    building and then coarsening. Holds for *any* fibre map, contiguous or not —
    no order-preservation and no surjectivity are needed.

    There is one side condition, `hNJ`: every coarse group must have nonzero
    total population. It is not implied by surjectivity, since a coarse group
    can be hit only by fine groups that are themselves empty. Under the
    nonnegative populations `ContactMatrix` enforces it *is* discharged, because a
    fibre summing to zero forces each of its members to zero — see
    `weightedCoarsen_proportionate_nonneg` immediately above, which is the version
    to cite for anything about ContACT's own domain.

    Distinct from `coarsenAction_comp`, which is functoriality of coarsening
    alone and holds for *any* matrix. This theorem relates two different
    operations — `build` and `coarsen` — and therefore depends on the specific
    lift used. It fails for the naive `a_i*a_j/D` form, which is not valid
    MeanContacts when populations differ, and it also fails for the assortative/
    disassortative Britton-Ball kernels (`src/activity_refinement.jl`) unless
    the fibre map is order-preserving — see the Julia regression test in
    `test/runtests.jl` for both counterexamples.

    Consequence: matrices generated from a parameter vector stay consistent
    across resolutions, so parameter-level edits (e.g. an intervention scaling
    one group's activity) can be defined at one partition and interpreted at a
    coarser one. -/
theorem weightedCoarsen_proportionate {m n : ℕ}
    (a N : Fin m → ℝ) (f : Fin m → Fin n)
    (hNJ : ∀ J : Fin n, coarsenPop N f J ≠ 0) :
    weightedCoarsen (proportionateMatrix a N) N f =
    proportionateMatrix (coarsenActivity a N f) (coarsenPop N f) := by
  have hD := totalActivity_coarsen a N f hNJ
  ext I J
  have hfun : (fun i j => proportionateMatrix a N i j * N j)
      = (fun i j => (a i * N i) * (a j * N j) / totalActivity a N) := by
    funext i j
    unfold proportionateMatrix
    ring
  show coarsenTotal (fun i j => proportionateMatrix a N i j * N j) f I J / coarsenPop N f J
     = proportionateMatrix (coarsenActivity a N f) (coarsenPop N f) I J
  rw [hfun, coarsenTotal_outer_div]
  show coarsenActivityMass a N f I * coarsenActivityMass a N f J / totalActivity a N
       / coarsenPop N f J
     = coarsenActivity a N f I * coarsenPop N f I * coarsenActivity a N f J
       / totalActivity (coarsenActivity a N f) (coarsenPop N f)
  rw [hD, coarsenActivity_apply, coarsenActivity_apply,
      div_mul_cancel₀ _ (hNJ I)]
  ring
