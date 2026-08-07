import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import ContACTProofs.Coarsening

/-!
*Source: `PopulationTransport.lean`*

# Population Transport Properties

Formalises `transport_population` (`src/population_transport.jl`), the
identity/composition-respecting counterpart to `reproject`
(`ContACTProofs.Reprojection`). `reproject` fails the composition law even
starting from a matrix already reciprocal at its own population
(`reproject_not_functorial`, `Reprojection.lean`), and fails the identity law
whenever the input is *not* already reciprocal at the target population — it is the
identity exactly on the matrices already reciprocal there, which is a subspace of
the population fibre, not the whole of it (`reproject_id_of_reciprocal`,
`Reprojection.lean`). `populationTransport` holds total contacts fixed and satisfies
both laws exactly, for every matrix, not just a reciprocal fibre. On strictly
positive populations it is right-multiplication by `diag(source)·diag(target)⁻¹`,
giving an action of the positive diagonal group; that description does not extend to
the boundary where Julia admits a zero source population, since the map is not
invertible there.

The population-transport formula:
    M'[i,j] = M[i,j] · source_j / target_j

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Identity: transporting onto the source population is the identity, for every M | ✅ |
| 2 | Composition is exact, for every M | ✅ |
| 3 | Total contacts preserved entrywise: M'[i,j]·target_j = M[i,j]·source_j | ✅ |
| 4 | Reciprocity is preserved when the source is already reciprocal at `source` | ✅ |
| 5 | Invertible: round trip recovers the original, for every M (no reciprocity needed) | ✅ |
| 6 | Commutes with population-weighted coarsening, for any fibre map (`_nonneg`: no coarse-population condition) | ✅ |
| 7 | The nonzero-target hypothesis in (6) is load-bearing (counterexample) | ✅ |

-/

/-- Transport a contact matrix from `source` to `target`, holding total contacts
`M[i,j]·source_j` fixed. -/
noncomputable def populationTransport {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ)
    (source target : Fin n → ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => M i j * source j / target j

/-- **Identity**: transporting a matrix onto its own population changes nothing. -/
theorem populationTransport_id {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop : Fin n → ℝ) (hpop : ∀ j, pop j ≠ 0) :
    populationTransport M pop pop = M := by
  ext i j
  have hj := hpop j
  simp only [populationTransport]
  field_simp

/-- **Composition** is exact: routing a transport through an intermediate
population agrees with transporting directly, unlike `reproject`
(`reproject_not_functorial`, `Reprojection.lean`). -/
theorem populationTransport_comp {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (N N' N'' : Fin n → ℝ)
    (hN' : ∀ j, N' j ≠ 0)
    -- Not load-bearing: the equality holds even at `N'' j = 0` (both sides reduce
    -- to `0` under Lean's `x / 0 = 0`). Kept to match `_transport_at`'s domain,
    -- which requires every target entry strictly positive.
    (hN'' : ∀ j, N'' j ≠ 0) :
    populationTransport (populationTransport M N N') N' N'' =
    populationTransport M N N'' := by
  ext i j
  have hj' := hN' j
  simp only [populationTransport]
  field_simp

/-- **Total contacts are preserved entrywise**: `M'[i,j]·target_j = M[i,j]·source_j`.
Stronger than `reproject`'s total-contacts result (`reproject_total_contacts`,
`Reprojection.lean`), which only holds after summing over all pairs. -/
theorem populationTransport_total_contacts {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (source target : Fin n → ℝ)
    (htarget : ∀ j, target j ≠ 0) (i j : Fin n) :
    populationTransport M source target i j * target j = M i j * source j := by
  have hj := htarget j
  simp only [populationTransport]
  field_simp

/-- **Reciprocity is preserved conditionally**: if `M` is already reciprocal at
`source`, its transport is reciprocal at `target`. Unlike `reproject`, which
repairs reciprocity unconditionally, transport only carries reciprocity forward —
it does not create it (see the accompanying Julia regression test for a
non-reciprocal counterexample where transport leaves the imbalance in place). -/
theorem populationTransport_reciprocity {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (source target : Fin n → ℝ)
    (htarget : ∀ j, target j ≠ 0)
    (hrec : ∀ i j, M i j * source j = M j i * source i)
    (i j : Fin n) :
    populationTransport M source target i j * target j =
    populationTransport M source target j i * target i := by
  rw [populationTransport_total_contacts M source target htarget i j,
      populationTransport_total_contacts M source target htarget j i]
  exact hrec i j

/-- **Invertible**: transporting to `N'` and back to `N` recovers the original
matrix exactly, for *every* `M` — no reciprocity hypothesis needed. (Reciprocity
only matters for `populationTransport_reciprocity` above, which is about what the
transported matrix satisfies, not about invertibility.) A corollary of
`populationTransport_comp` and `populationTransport_id` alone. -/
theorem populationTransport_roundtrip {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (N N' : Fin n → ℝ)
    (hN : ∀ j, N j ≠ 0) (hN' : ∀ j, N' j ≠ 0) :
    populationTransport (populationTransport M N N') N' N = M := by
  rw [populationTransport_comp M N N' N hN' hN, populationTransport_id M N hN]

/-- **Commutes with coarsening.** In total-contacts terms coarsening is a plain sum
over each fibre and `populationTransport` leaves total contacts alone, so the two
commute whenever the coarse target population is the pushforward `coarsenPop N' f`
of the fine one.

Note what is *not* assumed: nothing at all about `f`. The fibres need not be
contiguous or order-preserving, unlike the Britton–Ball assortative and
disassortative activity kernels, which allocate contacts by walking strata in order
and so commute with coarsening only for order-preserving fibre maps. The two
hypotheses are exactly the divisions that have to cancel — `N' j` on the left, where
coarsening multiplies back what transport divided out, and `coarsenPop N f J` on the
right, where transport multiplies back what coarsening divided out.

Models the Julia regression test `"transport: coarsening equivariance"`
(`test/runtests.jl`), which exercises non-contiguous fibre maps such as `[1,2,1,2]`
precisely because this theorem needs no order condition. -/
theorem weightedCoarsen_populationTransport {m n : ℕ}
    (M : Matrix (Fin m) (Fin m) ℝ) (N N' : Fin m → ℝ) (f : Fin m → Fin n)
    (hN' : ∀ j, N' j ≠ 0)
    (hNJ : ∀ J, coarsenPop N f J ≠ 0) :
    weightedCoarsen (populationTransport M N N') N' f =
    populationTransport (weightedCoarsen M N f) (coarsenPop N f) (coarsenPop N' f) := by
  ext I J
  have hcancel : ∀ i j, M i j * N j / N' j * N' j = M i j * N j :=
    fun i j => div_mul_cancel₀ _ (hN' j)
  simp only [weightedCoarsen, populationTransport, hcancel]
  rw [div_mul_cancel₀ _ (hNJ J)]

/-- The nonzero-target hypothesis in `weightedCoarsen_populationTransport` is
load-bearing, not defensive bookkeeping. Merging `{0, 1}` into a single coarse group
with `N = (1, 1)`, `N' = (0, 1)` and `M ≡ 1`, the fibre still has a survivor, so
`coarsenPop N' f` is nonzero and the other hypothesis holds. Transport zeroes the
`N' = 0` column before coarsening can sum it, so the left side loses that column's
contribution while the right side keeps it.

Contrast `hNJ`, which is load-bearing only because `N : Fin m → ℝ` admits negative
entries here: for the nonnegative populations `ContactMatrix` enforces,
`coarsenPop N f J = 0` forces every `N j` in the fibre to vanish, and both sides are
then `0`. -/
theorem weightedCoarsen_populationTransport_needs_nonzero_target :
    ∃ (M : Matrix (Fin 2) (Fin 2) ℝ) (N N' : Fin 2 → ℝ) (f : Fin 2 → Fin 1),
      (∀ J, coarsenPop N f J ≠ 0) ∧
      weightedCoarsen (populationTransport M N N') N' f ≠
        populationTransport (weightedCoarsen M N f) (coarsenPop N f) (coarsenPop N' f) := by
  refine ⟨fun _ _ => 1, fun _ => 1, ![0, 1], fun _ => 0, ?_, ?_⟩
  · intro J
    obtain rfl : J = 0 := Fin.fin_one_eq_zero J
    norm_num [coarsenPop, Fin.sum_univ_two]
  · intro h
    have := congrFun (congrFun h 0) 0
    simp [weightedCoarsen, populationTransport, coarsenPop, coarsenTotal] at this

/-- **Transport–coarsening equivariance on the domain the package actually uses.**
Identical conclusion to `weightedCoarsen_populationTransport`, with the
coarse-population hypothesis `hNJ` replaced by nonnegativity of the fine source
population — which every `ContactMatrix` satisfies by construction. A fibre summing
to zero then forces each of its members to zero, so the coarse total contacts vanish
too and both sides are `0`; that is exactly `weightedCoarsen_mul_pop`.

This is the theorem the README row and `transport_population`'s docstring describe
when they say the law holds for any fibre map. -/
theorem weightedCoarsen_populationTransport_nonneg {m n : ℕ}
    (M : Matrix (Fin m) (Fin m) ℝ) (N N' : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) (hN' : ∀ j, N' j ≠ 0) :
    weightedCoarsen (populationTransport M N N') N' f =
    populationTransport (weightedCoarsen M N f) (coarsenPop N f) (coarsenPop N' f) := by
  ext I J
  have hfun : (fun i j => populationTransport M N N' i j * N' j)
      = (fun i j => M i j * N j) := by
    funext i j
    simp only [populationTransport]
    exact div_mul_cancel₀ _ (hN' j)
  show coarsenTotal (fun i j => populationTransport M N N' i j * N' j) f I J
         / coarsenPop N' f J
     = weightedCoarsen M N f I J * coarsenPop N f J / coarsenPop N' f J
  rw [hfun, weightedCoarsen_mul_pop M N f hN I J]
