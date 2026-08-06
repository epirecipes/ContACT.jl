import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp

/-!
*Source: `PopulationTransport.lean`*

# Population Transport Properties

Formalises `transport_population` (`src/population_transport.jl`), the
identity/composition-respecting counterpart to `reproject`
(`ContACTProofs.Reprojection`). `reproject` fails the composition law even
starting from a matrix already reciprocal at its own population
(`reproject_not_functorial`, `Reprojection.lean`), and fails the identity law
whenever the input is *not* already reciprocal at the target population — it is
the identity only on that fibre (`reproject_id_of_reciprocal`, same file).
`populationTransport` holds total contacts fixed and satisfies both laws exactly,
for every matrix, not just a reciprocal fibre — it is right-multiplication by
`diag(source)·diag(target)⁻¹`, a diagonal-group action.

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
