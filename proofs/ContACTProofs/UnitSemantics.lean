import Mathlib.Data.Matrix.Basic
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.Data.Real.Basic
import ContACTProofs.Symmetrisation

/-!
*Source: `UnitSemantics.lean`*

# Naturality of unit semantics

`ContactCounts`, `MeanContacts` and `PerCapitaRate` are not three different objects:
they are three *representations* of one object, related by invertible rescalings of
the rows and columns. The quantity they all scale is total contacts
`T[i,j] = M[i,j] · N_j`, and a representation is a pair of exponents `(a, b)`:

    R[i,j] = T[i,j] · N_i^(-a) · N_j^(-b)

| Semantics       | `(a,b)` | entry          | row `r_i`  | col `c_j` |
|-----------------|---------|----------------|------------|-----------|
| `ContactCounts` | `(0,0)` | `M[i,j]·N_j`   | `1`        | `N_j`     |
| `MeanContacts`  | `(0,1)` | `M[i,j]`       | `1`        | `1`       |
| `PerCapitaRate` | `(1,1)` | `M[i,j]/N_i`   | `1 / N_i`  | `1`       |

writing `R = rescale M r c` relative to the `MeanContacts` matrix `M`.

A morphism applied directly in the wrong representation is not natural —
`symmetrise_not_natural_naive` is a concrete `2 × 2` counterexample. `transport`
conjugates a morphism written in its home representation, and `transport_naturality`
holds for *any* operation, which is how `src/unit_semantics.jl` makes `↔`, `↓`, `↑`
and `⊗` natural at once.

For symmetrisation the transport has a closed form: `symmetriseAt`, the ordinary
`MeanContacts` symmetrisation at the effective population `v_k = N_k · r_k / c_k`.
Every law of `symmetriseMatrix` therefore carries over, and `v` collapses to `1` for
both `ContactCounts` and `PerCapitaRate`, where symmetrisation is plain averaging.

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Row/column rescalings form a groupoid | ✅ |
| 2 | A morphism applied directly in the wrong representation is not natural | ✅ |
| 3 | `transport` is natural for any operation | ✅ |
| 4 | `transport` is the identity in the home representation | ✅ |
| 5 | `symmetriseAt` is the transport of `symmetriseMatrix` | ✅ |
| 6 | Closed forms under `ContactCounts` and `PerCapitaRate` (plain averaging) | ✅ |
| 7 | `S_Φ` establishes reciprocity, and is idempotent, in every representation | ✅ |

-/

-- ═══════════════════════════════════════════════════════════════════════════
-- The groupoid of representations
-- ═══════════════════════════════════════════════════════════════════════════

/-- Change of representation: multiply row `i` by `row i` and column `j` by
`col j`. This is the shape of every `UnitSemantics` conversion, which in the Julia
code all route through `reinterpret_units`. -/
def rescale {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (row col : Fin n → ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun i j => M i j * row i * col j

theorem rescale_apply {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (row col : Fin n → ℝ)
    (i j : Fin n) : rescale M row col i j = M i j * row i * col j := rfl

/-- Identity: converting to the representation a matrix is already in does
nothing. -/
theorem rescale_one {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) :
    rescale M (fun _ => 1) (fun _ => 1) = M := by
  ext i j; simp [rescale]

/-- Composition: representation changes compose by multiplying the row and column
scalings independently (equivalently, by adding the exponents `a` and `b`). -/
theorem rescale_rescale {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ)
    (row col row' col' : Fin n → ℝ) :
    rescale (rescale M row col) row' col' =
      rescale M (fun i => row i * row' i) (fun j => col j * col' j) := by
  ext i j; simp [rescale]; ring

/-- Invertibility: every representation change is an isomorphism, provided no
scaling vanishes. This is the domain restriction `reinterpret_units` enforces by
throwing on a zero-population row or column with nonzero contacts. -/
theorem rescale_inv {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (row col : Fin n → ℝ)
    (hrow : ∀ i, row i ≠ 0) (hcol : ∀ j, col j ≠ 0) :
    rescale (rescale M row col) (fun i => (row i)⁻¹) (fun j => (col j)⁻¹) = M := by
  ext i j
  have hri : row i ≠ 0 := hrow i
  have hcj : col j ≠ 0 := hcol j
  simp only [rescale]
  field_simp

-- ═══════════════════════════════════════════════════════════════════════════
-- Transport
-- ═══════════════════════════════════════════════════════════════════════════

/-- Conjugate an operation by a change of representation: undo the rescaling, apply
`f`, redo it. This is `_via` in `src/unit_semantics.jl`. -/
noncomputable def transport {n : ℕ} (row col : Fin n → ℝ)
    (f : Matrix (Fin n) (Fin n) ℝ → Matrix (Fin n) (Fin n) ℝ) :
    Matrix (Fin n) (Fin n) ℝ → Matrix (Fin n) (Fin n) ℝ :=
  fun R => rescale (f (rescale R (fun i => (row i)⁻¹) (fun j => (col j)⁻¹))) row col

/-- **Naturality, for any operation whatever**: `Φ ∘ f = transport(f) ∘ Φ`. No
hypothesis on `f` is needed, only that the rescalings are invertible. -/
theorem transport_naturality {n : ℕ} (row col : Fin n → ℝ)
    (hrow : ∀ i, row i ≠ 0) (hcol : ∀ j, col j ≠ 0)
    (f : Matrix (Fin n) (Fin n) ℝ → Matrix (Fin n) (Fin n) ℝ)
    (M : Matrix (Fin n) (Fin n) ℝ) :
    rescale (f M) row col = transport row col f (rescale M row col) := by
  rw [transport, rescale_inv M row col hrow hcol]

/-- Undoing a rescaling before applying it. -/
theorem rescale_inv' {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (row col : Fin n → ℝ)
    (hrow : ∀ i, row i ≠ 0) (hcol : ∀ j, col j ≠ 0) :
    rescale (rescale M (fun i => (row i)⁻¹) (fun j => (col j)⁻¹)) row col = M := by
  ext i j
  have hri : row i ≠ 0 := hrow i
  have hcj : col j ≠ 0 := hcol j
  simp only [rescale]
  field_simp

/-- In the home representation transport is the identity. -/
theorem transport_one {n : ℕ}
    (f : Matrix (Fin n) (Fin n) ℝ → Matrix (Fin n) (Fin n) ℝ)
    (M : Matrix (Fin n) (Fin n) ℝ) :
    transport (fun _ => 1) (fun _ => 1) f M = f M := by
  simp only [transport, inv_one, rescale_one]

-- ═══════════════════════════════════════════════════════════════════════════
-- Symmetrisation, transported to an arbitrary representation
-- ═══════════════════════════════════════════════════════════════════════════

/-!
If `R = rescale M row col` then total contacts satisfy
`T[i,j] = R[i,j] · (N_j · row_j / col_j) · (row_i / row_j)⁻¹`… more usefully, the
whole representation behaves like the `MeanContacts` one with population vector
`v_k = N_k · row_k / col_k`. That is the definition below, and it is what lets the
reciprocity and idempotence theorems of `Symmetrisation.lean` be reused verbatim.
-/

/-- Effective population of a representation: the vector that `symmetriseMatrix`
must be evaluated at for symmetrisation in this representation. -/
noncomputable def effectivePop {n : ℕ} (row col pop : Fin n → ℝ) : Fin n → ℝ :=
  fun k => pop k * row k / col k

/-- Symmetrisation in the representation given by row scaling `row` and column
scaling `col`. -/
noncomputable def symmetriseAt {n : ℕ} (row col pop : Fin n → ℝ)
    (R : Matrix (Fin n) (Fin n) ℝ) : Matrix (Fin n) (Fin n) ℝ :=
  symmetriseMatrix R (effectivePop row col pop)

theorem symmetriseAt_apply {n : ℕ} (row col pop : Fin n → ℝ)
    (R : Matrix (Fin n) (Fin n) ℝ) (i j : Fin n) :
    symmetriseAt row col pop R i j =
      (R i j * effectivePop row col pop j + R j i * effectivePop row col pop i) /
        (2 * effectivePop row col pop j) := rfl

/-- `MeanContacts` is `row = col = 1`, where this is the symmetrisation already
formalised in `Symmetrisation.lean`. -/
theorem symmetriseAt_meanContacts {n : ℕ} (pop : Fin n → ℝ)
    (R : Matrix (Fin n) (Fin n) ℝ) :
    symmetriseAt (fun _ => 1) (fun _ => 1) pop R = symmetriseMatrix R pop := by
  ext i j; simp [symmetriseAt, effectivePop, symmetriseMatrix]

-- ═══════════════════════════════════════════════════════════════════════════
-- 1. The naive operation is not natural
-- ═══════════════════════════════════════════════════════════════════════════

/-- **The latent bug, as a counterexample.** Applying the `MeanContacts` formula to
a matrix that is *not* in the `MeanContacts` representation does not commute with
the change of representation. Here `col = N` (the `to_counts` isomorphism), so the
left-hand side is what the old `symmetrise` returned when handed a `ContactCounts`
matrix, and the right-hand side is what naturality demands.

Witnesses: `N = (1, 2)` and the mean-contacts matrix with a single nonzero entry
`M[0,1] = 1`. The two sides disagree at entry `(1,0)` — `2` against `1`. -/
theorem symmetrise_not_natural_naive :
    ∃ (M : Matrix (Fin 2) (Fin 2) ℝ) (pop : Fin 2 → ℝ),
      (∀ i, 0 < pop i) ∧
      symmetriseMatrix (rescale M (fun _ => 1) pop) pop ≠
        rescale (symmetriseMatrix M pop) (fun _ => 1) pop := by
  refine ⟨!![0, 1; 0, 0], ![1, 2], ?_, ?_⟩
  · intro i; fin_cases i <;> norm_num
  · intro h
    have h10 := congrFun (congrFun h 1) 0
    simp [symmetriseMatrix, rescale] at h10

-- ═══════════════════════════════════════════════════════════════════════════
-- 2. Naturality of the representation-aware operator
-- ═══════════════════════════════════════════════════════════════════════════

/-- **Naturality square.** `Φ ∘ S = S_Φ ∘ Φ`, where `Φ = rescale · row col` is the
change of representation, `S = symmetriseMatrix` is symmetrisation in
`MeanContacts`, and `S_Φ = symmetriseAt row col` is symmetrisation in the target
representation.

In Julia this is
`reinterpret_units(symmetrise(cm), s) == symmetrise(reinterpret_units(cm, s))`,
which holds because `symmetrise` dispatches on `cm.semantics`. -/
theorem symmetrise_naturality {n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (pop row col : Fin n → ℝ)
    (hpop : ∀ i, pop i ≠ 0) (hrow : ∀ i, row i ≠ 0) (hcol : ∀ j, col j ≠ 0) :
    rescale (symmetriseMatrix M pop) row col =
      symmetriseAt row col pop (rescale M row col) := by
  ext i j
  have hpi : pop i ≠ 0 := hpop i
  have hpj : pop j ≠ 0 := hpop j
  have hri : row i ≠ 0 := hrow i
  have hrj : row j ≠ 0 := hrow j
  have hci : col i ≠ 0 := hcol i
  have hcj : col j ≠ 0 := hcol j
  rw [rescale_apply, symmetriseAt_apply, rescale_apply, rescale_apply]
  simp only [symmetriseMatrix, effectivePop]
  field_simp

/-- `symmetriseAt` is the transport of `symmetriseMatrix`: the closed form and the
conjugation define the same operation. -/
theorem symmetriseAt_eq_transport {n : ℕ} (row col pop : Fin n → ℝ)
    (hpop : ∀ i, pop i ≠ 0) (hrow : ∀ i, row i ≠ 0) (hcol : ∀ j, col j ≠ 0)
    (R : Matrix (Fin n) (Fin n) ℝ) :
    symmetriseAt row col pop R =
      transport row col (fun M => symmetriseMatrix M pop) R := by
  rw [transport, symmetrise_naturality _ pop row col hpop hrow hcol,
      rescale_inv' R row col hrow hcol]

-- ═══════════════════════════════════════════════════════════════════════════
-- 3. Closed form in each representation
-- ═══════════════════════════════════════════════════════════════════════════

/-- Whenever the effective population is `1`, symmetrisation is plain averaging.
This is the `d = a - b = 0` case: the reciprocity law reduces to plain matrix
symmetry, so no population survives in the formula. -/
theorem symmetriseAt_of_effectivePop_one {n : ℕ} (row col pop : Fin n → ℝ)
    (R : Matrix (Fin n) (Fin n) ℝ) (h : ∀ k, effectivePop row col pop k = 1)
    (i j : Fin n) :
    symmetriseAt row col pop R i j = (R i j + R j i) / 2 := by
  rw [symmetriseAt_apply, h i, h j]; ring

/-- **`ContactCounts`** is `row = 1`, `col = N`: the effective population is `1`, so
symmetrisation of total contacts is plain symmetrisation of the matrix. -/
theorem symmetriseAt_counts {n : ℕ} (R : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ) (hpop : ∀ i, pop i ≠ 0) (i j : Fin n) :
    symmetriseAt (fun _ => 1) pop pop R i j = (R i j + R j i) / 2 := by
  refine symmetriseAt_of_effectivePop_one _ _ _ _ (fun k => ?_) i j
  simp only [effectivePop, mul_one]
  exact div_self (hpop k)

/-- **`PerCapitaRate`** is `row = 1/N`, `col = 1`: the effective population is again
`1`. This is the representation in which a reciprocal matrix is *symmetric* — the
`β_ij` of Wallinga et al. (2006) and Lomas et al. (2025), and `socialmixr`'s
`per_capita`. -/
theorem symmetriseAt_perCapita {n : ℕ} (R : Matrix (Fin n) (Fin n) ℝ)
    (pop : Fin n → ℝ) (hpop : ∀ i, pop i ≠ 0) (i j : Fin n) :
    symmetriseAt (fun k => (pop k)⁻¹) (fun _ => 1) pop R i j = (R i j + R j i) / 2 := by
  refine symmetriseAt_of_effectivePop_one _ _ _ _ (fun k => ?_) i j
  simp only [effectivePop, div_one]
  exact mul_inv_cancel₀ (hpop k)

-- ═══════════════════════════════════════════════════════════════════════════
-- 4. The laws survive transport
-- ═══════════════════════════════════════════════════════════════════════════

/-- Reciprocity in a representation: total contacts are symmetric, expressed
through the effective population. At `row = col = 1` this is
`M[i,j]·N_j = M[j,i]·N_i`; for `ContactCounts` and `PerCapitaRate` the effective
population is `1` and it degenerates to plain symmetry of the matrix. -/
def ReciprocalAt {n : ℕ} (row col pop : Fin n → ℝ)
    (R : Matrix (Fin n) (Fin n) ℝ) : Prop :=
  ∀ i j, R i j * effectivePop row col pop j = R j i * effectivePop row col pop i

/-- `symmetriseAt` establishes reciprocity in its own representation. -/
theorem symmetriseAt_reciprocity {n : ℕ} (R : Matrix (Fin n) (Fin n) ℝ)
    (row col pop : Fin n → ℝ) (h : ∀ i, effectivePop row col pop i > 0) :
    ReciprocalAt row col pop (symmetriseAt row col pop R) :=
  fun i j => symmetrise_reciprocity R (effectivePop row col pop) h i j

/-- `symmetriseAt` is idempotent in every representation. -/
theorem symmetriseAt_idempotent {n : ℕ} (R : Matrix (Fin n) (Fin n) ℝ)
    (row col pop : Fin n → ℝ) (h : ∀ i, effectivePop row col pop i > 0) :
    symmetriseAt row col pop (symmetriseAt row col pop R) = symmetriseAt row col pop R :=
  symmetrise_idempotent R (effectivePop row col pop) h

/-- In a representation whose effective population is `1` — `ContactCounts` and
`PerCapitaRate` — reciprocity *is* symmetry of the matrix. This is the property
`PerCapitaRate` is named for. -/
theorem reciprocalAt_iff_symm {n : ℕ} (row col pop : Fin n → ℝ)
    (R : Matrix (Fin n) (Fin n) ℝ) (h : ∀ k, effectivePop row col pop k = 1) :
    ReciprocalAt row col pop R ↔ ∀ i j, R i j = R j i := by
  constructor <;> intro hR i j
  · have := hR i j; rw [h i, h j] at this; simpa using this
  · rw [h i, h j]; simpa using hR i j
