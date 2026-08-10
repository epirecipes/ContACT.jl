import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import ContACTProofs.Coarsening

/-!
*Source: `Refinement.lean`*

# Refinement as a Section of Coarsening

Formalises `refine` (`src/refinement.jl`), which splits each coarse group's contacts
between its fine subgroups in proportion to their populations:

    M_fine[i,j] = M_coarse[f i, f j] · N_i / N_{f i}

Only the contactee index carries a population share. In total contacts the split is
proportional on both axes; the participant-side share cancels against the division by
participant population that turns total contacts into mean contacts.

The central result is that refinement is a **section** of the population-weighted
coarsening of `ContACTProofs.Coarsening`: coarsening a refined matrix returns the
original. It is not a retraction — refining a coarsened matrix does not recover it,
since coarsening discards the within-group structure that the prior then has to
guess (`refine_not_left_inverse`).

An empty coarse group supplies no shares to split by, so the section law needs its
row and column to be zero — the condition `refine` enforces at construction. That is
stricter than the neighbouring zero-population guards, not the same as them:
`transport_population` checks the participant column only, and `reinterpret_units`
checks whichever axis it is rescaling, while refinement needs both. Where no coarse
population vanishes there is no side condition at all
(`weightedCoarsen_refine_pos`).

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Section law: coarsening a refined matrix returns the original | ✅ |
| 2 | Section law with no side condition where no coarse population is zero | ✅ |
| 3 | Total contacts are preserved | ✅ |
| 4 | Reciprocity at the coarse population is inherited by the refined matrix | ✅ |
| 5 | Refinement is not a left inverse of coarsening (counterexample) | ✅ |

-/

/-- The share of its coarse group's contacts that fine group `i` receives,
    `N_i / N_{f i}`. An empty coarse group gives `0 / 0 = 0`, which is the zero
    branch `refine` takes in `src/refinement.jl`. -/
noncomputable def refineShare {m n : ℕ}
    (N : Fin m → ℝ) (f : Fin m → Fin n) (i : Fin m) : ℝ :=
  N i / coarsenPop N f (f i)

/-- Population-proportional refinement, matching `refine` in `src/refinement.jl`.
    A fine group with zero population gets a zero participant column. -/
noncomputable def refineMatrix {m n : ℕ}
    (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n) :
    Matrix (Fin m) (Fin m) ℝ :=
  fun i j => if N j = 0 then 0 else Mc (f i) (f j) * refineShare N f i

/-- The zero-column branch is invisible in total contacts: it fires only where the
    participant population, and hence the entry's contribution, is zero. Every result
    below goes through this, so none of them depends on that branch. -/
theorem refineMatrix_mul_pop {m n : ℕ}
    (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n) (i j : Fin m) :
    refineMatrix Mc N f i j * N j = Mc (f i) (f j) * refineShare N f i * N j := by
  unfold refineMatrix
  split_ifs with h
  · rw [h]; ring
  · rfl

/-- The shares within a coarse group sum to one, except at an empty group, where
    every share is `0 / 0 = 0` and they sum to zero. -/
theorem sum_refineShare {m n : ℕ}
    (N : Fin m → ℝ) (f : Fin m → Fin n) (I : Fin n) :
    (∑ i ∈ Finset.univ.filter (fun i => f i = I), refineShare N f i) =
    coarsenPop N f I / coarsenPop N f I := by
  have h : ∀ i ∈ Finset.univ.filter (fun i => f i = I),
      refineShare N f i = N i / coarsenPop N f I := by
    intro i hi
    rw [Finset.mem_filter] at hi
    unfold refineShare
    rw [hi.2]
  rw [Finset.sum_congr rfl h]
  simp only [div_eq_mul_inv]
  rw [← Finset.sum_mul]
  rfl

/-- **Refinement is a section of coarsening**: coarsening a refined matrix along the
    same fibre map returns the coarse matrix it came from.

    The two hypotheses are the guard `refine` applies: a coarse group of zero
    population must carry no contacts, in its row (there are no shares to split them
    by) or its column (its fine subgroups are empty, so they cannot carry it). Nothing
    is assumed about the fibre map — not surjectivity, and not order-preservation. -/
theorem weightedCoarsen_refine {m n : ℕ}
    (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hrow : ∀ I J, coarsenPop N f I = 0 → Mc I J = 0)
    (hcol : ∀ I J, coarsenPop N f J = 0 → Mc I J = 0) :
    weightedCoarsen (refineMatrix Mc N f) N f = Mc := by
  ext I J
  rw [weightedCoarsen_apply]
  have hblock : coarsenTotal (fun i j => refineMatrix Mc N f i j * N j) f I J
      = (∑ i ∈ Finset.univ.filter (fun i => f i = I), refineShare N f i)
        * (Mc I J * coarsenPop N f J) := by
    rw [Finset.sum_mul]
    simp only [coarsenTotal]
    refine Finset.sum_congr rfl fun i hi => ?_
    rw [Finset.mem_filter] at hi
    rw [coarsenPop, Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [Finset.mem_filter] at hj
    rw [refineMatrix_mul_pop, hi.2, hj.2]
    ring
  rw [hblock, sum_refineShare]
  rcases eq_or_ne (coarsenPop N f J) 0 with hJ | hJ
  · rw [hJ, mul_zero, mul_zero, zero_div, hcol I J hJ]
  · rcases eq_or_ne (coarsenPop N f I) 0 with hI | hI
    · rw [hI, div_zero, zero_mul, zero_div, hrow I J hI]
    · rw [div_self hI, one_mul, mul_div_assoc, div_self hJ, mul_one]

/-- The section law where no coarse population vanishes, which makes the empty-group
    hypotheses vacuous. The hypothesis is `≠ 0` rather than `0 <`; since
    `ContactMatrix` rejects negative populations, the two coincide on every matrix the
    package can build. This is the form to cite for any matrix whose groups are all
    populated. -/
theorem weightedCoarsen_refine_pos {m n : ℕ}
    (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hNP : ∀ I, coarsenPop N f I ≠ 0) :
    weightedCoarsen (refineMatrix Mc N f) N f = Mc :=
  weightedCoarsen_refine Mc N f
    (fun I _ h => absurd h (hNP I)) (fun _ J h => absurd h (hNP J))

/-- **Total contacts are preserved** by refinement: the fine total contacts equal the
    coarse ones. A corollary of the section law and
    `weightedCoarsen_preserves_total`. -/
theorem refine_preserves_total {m n : ℕ}
    (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j)
    (hrow : ∀ I J, coarsenPop N f I = 0 → Mc I J = 0)
    (hcol : ∀ I J, coarsenPop N f J = 0 → Mc I J = 0) :
    (∑ i : Fin m, ∑ j : Fin m, refineMatrix Mc N f i j * N j) =
    ∑ I : Fin n, ∑ J : Fin n, Mc I J * coarsenPop N f J := by
  rw [← weightedCoarsen_preserves_total (refineMatrix Mc N f) N f hN,
      weightedCoarsen_refine Mc N f hrow hcol]

/-- Under coarse reciprocity, a coarse entry divided by its *contactee* population is
    symmetric. This is the coarse matrix in `PerCapitaRate`, and it is what makes the
    refined matrix reciprocal: the refinement multiplies it by `N_i · N_j`, which is
    symmetric in `i` and `j`. -/
theorem coarse_ratio_symm {m n : ℕ}
    (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hcol : ∀ I J, coarsenPop N f J = 0 → Mc I J = 0)
    (hrec : ∀ I J, Mc I J * coarsenPop N f J = Mc J I * coarsenPop N f I)
    (I J : Fin n) :
    Mc I J / coarsenPop N f I = Mc J I / coarsenPop N f J := by
  rcases eq_or_ne (coarsenPop N f I) 0 with hI | hI
  · rw [hI, div_zero, hcol J I hI, zero_div]
  · rcases eq_or_ne (coarsenPop N f J) 0 with hJ | hJ
    · rw [hJ, div_zero, hcol I J hJ, zero_div]
    · field_simp
      linarith [hrec I J]

/-- **Reciprocity is inherited**: if the coarse matrix is reciprocal at the coarse
    population, the refined matrix is reciprocal at the fine one. Only the
    zero-column half of the `refine` guard is needed. -/
theorem refine_reciprocity {m n : ℕ}
    (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hcol : ∀ I J, coarsenPop N f J = 0 → Mc I J = 0)
    (hrec : ∀ I J, Mc I J * coarsenPop N f J = Mc J I * coarsenPop N f I)
    (i j : Fin m) :
    refineMatrix Mc N f i j * N j = refineMatrix Mc N f j i * N i := by
  rw [refineMatrix_mul_pop, refineMatrix_mul_pop]
  unfold refineShare
  have h := coarse_ratio_symm Mc N f hcol hrec (f i) (f j)
  calc Mc (f i) (f j) * (N i / coarsenPop N f (f i)) * N j
      = Mc (f i) (f j) / coarsenPop N f (f i) * (N i * N j) := by ring
    _ = Mc (f j) (f i) / coarsenPop N f (f j) * (N i * N j) := by rw [h]
    _ = Mc (f j) (f i) * (N j / coarsenPop N f (f j)) * N i := by ring

/-- Refinement is **not** a left inverse of coarsening. Two equally populated fine
    groups inside one coarse group, with a purely assortative fine matrix: coarsening
    averages the diagonal and off-diagonal away, and the proportional prior returns a
    uniform matrix. The coarse matrix alone no longer determines the within-group
    structure, so no population prior recovers it — which is why `↑` is documented as
    a section of `↓` and not an inverse. -/
theorem refine_not_left_inverse :
    ∃ (M : Matrix (Fin 2) (Fin 2) ℝ) (N : Fin 2 → ℝ) (f : Fin 2 → Fin 1),
      refineMatrix (weightedCoarsen M N f) N f ≠ M := by
  refine ⟨fun i j => if i = j then 1 else 0, fun _ => 1, fun _ => 0, fun h => ?_⟩
  have h00 := congrFun (congrFun h 0) 0
  simp [refineMatrix, refineShare, weightedCoarsen, coarsenTotal, coarsenPop] at h00
