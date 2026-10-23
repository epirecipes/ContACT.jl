import Mathlib.CategoryTheory.Category.Basic
import Mathlib.CategoryTheory.Functor.Basic
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.Tactic.FieldSimp
import Mathlib.Algebra.BigOperators.Field
import ContACTProofs.Coarsening
import ContACTProofs.Symmetrisation
import ContACTProofs.Refinement
import ContACTProofs.ConstrainedLift

/-!
*Source: `ContactCat.lean`*

# Contact Structures in Total-Contact Coordinates

A **contact structure** is a triple `(X, N, C)`: a finite set of groups `X = Fin n`,
their populations `N`, and the total contacts `C[i,j] = M[i,j]·N_j` between contactee
group `i` (row) and participant group `j` (column). A partition map `f : X → Y` is a
**morphism** `(X, N, C) → (Y, N', C')` exactly when `N'` and `C'` are the pushforwards
of `N` and `C` along `f`:

    N'[J]   = Σ_{f(j)=J} N_j
    C'[I,J] = Σ_{f(i)=I, f(j)=J} C[i,j]

These form a category. Nothing is assumed about the map or the data: the map need not
be surjective, and populations and contacts may be any real numbers. The category laws
are the functoriality of the two pushforwards (`coarsen_id`, `coarsen_functorial`,
`coarsenPop_id`, `coarsenPop_comp`). A morphism is determined by its map, and every
structure has exactly one morphism out of it along each map, landing on its pushforward
(`ContactHom.cod_eq_push`).

The shipped operations are stated against this category as follows.

- **`↓` (`coarsen`).** In `MeanContacts`, `weightedCoarsen` computes the target of the
  morphism along `f` (`ofMean_weightedCoarsen`), on nonnegative populations, including
  empty groups and non-surjective maps. A mean-contacts matrix represents its contact
  structure faithfully only in **normal form**, where every empty participant group has
  a zero column: coarsening along `id` is the projection onto normal form
  (`weightedCoarsen_id`), and is the identity exactly there. Normal form holds whenever
  every population is positive, and for every output of `coarsen`.
- **`↔` (`symmetrise`).** In total contacts, symmetrisation is `(C + Cᵀ)/2` at a fixed
  population. It is an endofunctor that leaves every map unchanged (`symFunctor`), so it
  commutes with pushforward (`ContactObj.push_sym`). In `MeanContacts` this gives
  `↓ ∘ ↔ = ↔ ∘ ↓`, empty groups included, wherever `symmetrise` accepts the input
  (`weightedCoarsen_symmetriseSafe`).
- **`⊕` (`compose_matrices`).** At a fixed population, total contacts form a commutative
  monoid under addition, and pushforward is a monoid homomorphism (`coarsenTotalHom`).
- **`↑` (`refine`).** Refinement constructs a finer structure together with a morphism
  back to the input matrix at the prior's aggregated population, which is the input
  itself when the prior aggregates exactly (`refineHom`).
- **`⊠` with a `GeneralizedLift` (`generalize`).** A block lift whose kernels sum to one
  on every block that carries contacts constructs a finer structure together with a
  morphism back to its base totals at the pushforward population (`blockLiftHom`).
  `generalize` lifts the symmetrised base totals, so in `MeanContacts` its coarsening is
  `symmetriseSafe` of the base (`weightedCoarsen_generalize`). The constrained lift
  (`⊠` with a `ConstrainedGeneralizedLift`) is not covered here.

Population-changing operations (`transport_population`, `reproject`) are not morphisms:
a morphism along `id` forces the population to stay the same.

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Contact structures and pushforward morphisms form a category, with no condition on the map or the data | ✅ |
| 2 | A morphism is determined by its map, and lands on the pushforward of its source | ✅ |
| 3 | Shipped `weightedCoarsen` computes the target of the morphism along a map (nonnegative populations; empty groups and non-surjective maps included) | ✅ |
| 4 | In `MeanContacts`, coarsening along `id` is the projection onto normal form; it is not the identity off normal form (counterexample) | ✅ |
| 5 | Symmetrisation is an endofunctor that leaves every map unchanged, commutes with pushforward, and is idempotent on objects | ✅ |
| 6 | `↓ ∘ ↔ = ↔ ∘ ↓` for the shipped operations, including empty groups, wherever `symmetrise` accepts the input; symmetrising after coarsening is defined on strictly more inputs | ✅ |
| 7 | Pushforward is an additive monoid homomorphism; weighted coarsening is additive with no hypothesis | ✅ |
| 8 | Refinement constructs a finer structure with a morphism back to the input matrix at the prior's aggregated population | ✅ |
| 9 | A normalised block lift constructs a finer structure with a morphism back to its base totals, and is reciprocal when its base totals are symmetric and each block's kernel is the transpose of its mirror block's | ✅ |
| 10 | In `MeanContacts`, coarsening a normalised block lift of symmetrised base totals returns `symmetriseSafe` of the base, which is the base when it is reciprocal and in normal form | ✅ |
| 11 | Column support and two-sided support are preserved by pushforward, and two-sided support by symmetrisation; on a mean-contacts structure, two-sided support is exactly the domain of `symmetrise` | ✅ |

-/

open CategoryTheory

-- ═══════════════════════════════════════════════════════════════════════════
-- The category
-- ═══════════════════════════════════════════════════════════════════════════

/-- Pushing populations forward along `id` leaves them unchanged. -/
theorem coarsenPop_id {n : ℕ} (N : Fin n → ℝ) : coarsenPop N id = N := by
  funext J
  simp only [coarsenPop, id]
  rw [show Finset.univ.filter (fun j => j = J) = {J} from by ext x; simp]
  simp

/-- A contact structure: groups `Fin n`, populations `pop`, and total contacts
    `counts` (row = contactee, column = participant). No sign or positivity
    condition is imposed. -/
structure ContactObj where
  n : ℕ
  pop : Fin n → ℝ
  counts : Matrix (Fin n) (Fin n) ℝ

/-- A morphism of contact structures: a partition map whose target data is the
    pushforward of its source data. -/
structure ContactHom (A B : ContactObj) where
  map : Fin A.n → Fin B.n
  pop_eq : coarsenPop A.pop map = B.pop
  counts_eq : coarsenTotal A.counts map = B.counts

/-- A morphism is determined by its map. -/
@[ext] theorem ContactHom.ext {A B : ContactObj} {φ ψ : ContactHom A B}
    (h : φ.map = ψ.map) : φ = ψ := by
  cases φ; cases ψ; cases h; rfl

instance : Category ContactObj where
  Hom := ContactHom
  id A := ⟨id, coarsenPop_id A.pop, coarsen_id A.counts⟩
  comp φ ψ :=
    ⟨ψ.map ∘ φ.map,
     by funext K; rw [← coarsenPop_comp, φ.pop_eq, ψ.pop_eq],
     by rw [coarsen_functorial, φ.counts_eq, ψ.counts_eq]⟩
  id_comp _ := ContactHom.ext rfl
  comp_id _ := ContactHom.ext rfl
  assoc _ _ _ := ContactHom.ext rfl

@[simp] theorem ContactHom.id_map (A : ContactObj) : (𝟙 A : A ⟶ A).map = id := rfl

@[simp] theorem ContactHom.comp_map {A B C : ContactObj} (φ : A ⟶ B) (ψ : B ⟶ C) :
    (φ ≫ ψ).map = ψ.map ∘ φ.map := rfl

/-- The pushforward of a contact structure along a map. -/
noncomputable def ContactObj.push (A : ContactObj) {m : ℕ} (f : Fin A.n → Fin m) :
    ContactObj :=
  ⟨m, coarsenPop A.pop f, coarsenTotal A.counts f⟩

/-- The morphism along `f` from a contact structure to its pushforward. -/
noncomputable def ContactObj.pushHom (A : ContactObj) {m : ℕ} (f : Fin A.n → Fin m) :
    A ⟶ A.push f :=
  ⟨f, rfl, rfl⟩

/-- Every morphism lands on the pushforward of its source along its map, so there is
    exactly one morphism out of a structure along each map. -/
theorem ContactHom.cod_eq_push {A B : ContactObj} (φ : A ⟶ B) : B = A.push φ.map := by
  obtain ⟨n, pop, counts⟩ := B
  obtain ⟨map, hp, hc⟩ := φ
  show ContactObj.mk n pop counts =
    ContactObj.mk n (coarsenPop A.pop map) (coarsenTotal A.counts map)
  rw [hp, show coarsenTotal A.counts map = counts from hc]

-- ═══════════════════════════════════════════════════════════════════════════
-- The MeanContacts representation and the shipped `coarsen`
-- ═══════════════════════════════════════════════════════════════════════════

/-- Mean contacts to total contacts: `C[i,j] = M[i,j]·N_j`. -/
noncomputable def toCounts {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun i j => M i j * N j

/-- Total contacts to mean contacts: `M[i,j] = C[i,j]/N_j`, which is `0` where
    `N_j = 0`. -/
noncomputable def ofCounts {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun i j => C i j / N j

/-- The contact structure of a mean-contacts matrix `M` with populations `N`. -/
noncomputable def ofMean (n : ℕ) (M : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ) :
    ContactObj :=
  ⟨n, N, toCounts M N⟩

/-- `weightedCoarsen` is, by definition, the mean-contacts form of the pushforward. -/
theorem weightedCoarsen_eq_ofCounts {m n : ℕ} (M : Matrix (Fin m) (Fin m) ℝ)
    (N : Fin m → ℝ) (f : Fin m → Fin n) :
    weightedCoarsen M N f = ofCounts (coarsenTotal (toCounts M N) f) (coarsenPop N f) :=
  rfl

/-- On nonnegative populations, `f` is a morphism from the structure of `(M, N)` to the
    structure of `(weightedCoarsen M N f, f_*N)`. Nothing is assumed about `f`, and
    empty groups are allowed. -/
noncomputable def coarsenHom {m n : ℕ} (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ)
    (f : Fin m → Fin n) (hN : ∀ j, 0 ≤ N j) :
    ofMean m M N ⟶ ofMean n (weightedCoarsen M N f) (coarsenPop N f) :=
  ⟨f, rfl, by
    funext I J
    exact (weightedCoarsen_mul_pop M N f hN I J).symm⟩

/-- **The shipped `coarsen` computes the pushforward**: on nonnegative populations, the
    structure of the coarsened matrix is the pushforward of the structure of the fine
    one. -/
theorem ofMean_weightedCoarsen {m n : ℕ} (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ)
    (f : Fin m → Fin n) (hN : ∀ j, 0 ≤ N j) :
    ofMean n (weightedCoarsen M N f) (coarsenPop N f) = (ofMean m M N).push f :=
  ContactHom.cod_eq_push (coarsenHom M N f hN)

/-- Converting to total contacts and back is the projection onto normal form: the
    column of an empty participant group becomes zero. -/
theorem ofCounts_toCounts {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ)
    (i j : Fin n) :
    ofCounts (toCounts M N) N i j = if N j = 0 then 0 else M i j := by
  unfold ofCounts toCounts
  split_ifs with h
  · rw [h, mul_zero, div_zero]
  · field_simp

/-- Converting to mean contacts and back loses nothing on total contacts that vanish at
    every empty participant group. -/
theorem toCounts_ofCounts {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ)
    (h : ∀ i j, N j = 0 → C i j = 0) : toCounts (ofCounts C N) N = C := by
  funext i j
  unfold toCounts ofCounts
  by_cases hj : N j = 0
  · rw [hj, mul_zero, h i j hj]
  · field_simp

/-- **The identity law in `MeanContacts`**: coarsening along `id` is the projection onto
    normal form. -/
theorem weightedCoarsen_id {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ) :
    weightedCoarsen M N id = fun i j => if N j = 0 then 0 else M i j := by
  funext i j
  rw [weightedCoarsen_eq_ofCounts, coarsen_id, coarsenPop_id, ofCounts_toCounts]

/-- Coarsening along `id` is the identity on matrices in normal form. -/
theorem weightedCoarsen_id_of_normal {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ)
    (N : Fin n → ℝ) (hnf : ∀ i j, N j = 0 → M i j = 0) :
    weightedCoarsen M N id = M := by
  funext i j
  rw [weightedCoarsen_id]
  dsimp only
  split_ifs with h
  · exact (hnf i j h).symm
  · rfl

/-- …and only there: a nonzero column at an empty group is zeroed. -/
theorem weightedCoarsen_id_not_id :
    ∃ (M : Matrix (Fin 2) (Fin 2) ℝ) (N : Fin 2 → ℝ),
      (∀ j, 0 ≤ N j) ∧ weightedCoarsen M N id ≠ M := by
  refine ⟨fun _ _ => 1, fun j => if j = 0 then 1 else 0, ?_, ?_⟩
  · intro j; by_cases h : j = 0 <;> simp [h]
  · intro h
    have h01 := congrFun (congrFun h 0) 1
    rw [weightedCoarsen_id] at h01
    simp at h01

-- ═══════════════════════════════════════════════════════════════════════════
-- `↔`: symmetrisation is an endofunctor over the maps
-- ═══════════════════════════════════════════════════════════════════════════

/-- Symmetrisation in total contacts: `(C + Cᵀ)/2`. -/
noncomputable def symCounts {n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ) :
    Matrix (Fin n) (Fin n) ℝ :=
  fun i j => (C i j + C j i) / 2

/-- Pushforward commutes with symmetrisation in total contacts, for any matrix. -/
theorem coarsenTotal_symCounts {m n : ℕ} (C : Matrix (Fin m) (Fin m) ℝ)
    (f : Fin m → Fin n) :
    coarsenTotal (symCounts C) f = symCounts (coarsenTotal C f) := by
  ext I J
  simp only [coarsenTotal, symCounts, add_div, Finset.sum_add_distrib, ← Finset.sum_div]
  congr 2
  exact Finset.sum_comm

/-- The symmetrisation of a contact structure, at the same population. -/
noncomputable def ContactObj.sym (A : ContactObj) : ContactObj :=
  ⟨A.n, A.pop, symCounts A.counts⟩

/-- A morphism of contact structures is also a morphism of their symmetrisations. -/
noncomputable def ContactHom.sym {A B : ContactObj} (φ : A ⟶ B) : A.sym ⟶ B.sym :=
  ⟨φ.map, φ.pop_eq, by
    show coarsenTotal (symCounts A.counts) φ.map = symCounts B.counts
    rw [coarsenTotal_symCounts, φ.counts_eq]⟩

/-- **Symmetrisation is an endofunctor** of contact structures. -/
noncomputable def symFunctor : ContactObj ⥤ ContactObj where
  obj := ContactObj.sym
  map := ContactHom.sym
  map_id _ := ContactHom.ext rfl
  map_comp _ _ := ContactHom.ext rfl

/-- …which leaves every map unchanged. -/
@[simp] theorem symFunctor_map_map {A B : ContactObj} (φ : A ⟶ B) :
    (symFunctor.map φ).map = φ.map := rfl

/-- Symmetrisation commutes with pushforward. -/
theorem ContactObj.push_sym (A : ContactObj) {m : ℕ} (f : Fin A.n → Fin m) :
    (A.push f).sym = A.sym.push f :=
  ContactHom.cod_eq_push (ContactHom.sym (A.pushHom f))

/-- Symmetrisation is idempotent on objects. -/
theorem ContactObj.sym_sym (A : ContactObj) : A.sym.sym = A.sym := by
  obtain ⟨n, pop, counts⟩ := A
  show ContactObj.mk n pop (symCounts (symCounts counts)) =
    ContactObj.mk n pop (symCounts counts)
  congr 1
  funext i j; simp only [symCounts]; ring

/-- **`↓ ∘ ↔ = ↔ ∘ ↓` for the shipped operations.** Coarsening commutes with
    symmetrisation in `MeanContacts`, for any map and including empty groups.

    `hcons` is the condition under which Julia's `symmetrise` accepts the fine
    matrix (`src/symmetrise.jl`, the zero-population branch of `_symmetrise_at`): no
    participant group reports contacts with an empty contactee group, so the empty
    group's *row* of total contacts is zero. Its column needs no condition. Under
    `hcons`, Julia's `symmetrise` also accepts the coarsened matrix
    (`weightedCoarsen_hcons`), and symmetrising after coarsening is defined on strictly
    more inputs (`weightedCoarsen_hcons_strict`). -/
theorem weightedCoarsen_symmetriseSafe {m n : ℕ}
    (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) (hcons : ∀ a b, N a = 0 → M a b * N b = 0) :
    weightedCoarsen (symmetriseSafe M N) N f =
      symmetriseSafe (weightedCoarsen M N f) (coarsenPop N f) := by
  funext I J
  have hcounts : toCounts (symmetriseSafe M N) N = symCounts (toCounts M N) := by
    funext i j
    exact symmetriseSafe_mul_pop M N hcons i j
  rw [weightedCoarsen_eq_ofCounts, hcounts, coarsenTotal_symCounts]
  unfold symmetriseSafe
  by_cases hJ : coarsenPop N f J = 0
  · rw [if_pos hJ]
    unfold ofCounts
    rw [hJ, div_zero]
  · rw [if_neg hJ, weightedCoarsen_mul_pop M N f hN I J, weightedCoarsen_mul_pop M N f hN J I]
    unfold ofCounts symCounts toCounts
    field_simp

/-- `hcons` passes to the coarsened matrix at the pushforward population: wherever
    Julia's `symmetrise` accepts a matrix, it also accepts that matrix's coarsening. -/
theorem weightedCoarsen_hcons {m n : ℕ} (M : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ)
    (f : Fin m → Fin n) (hN : ∀ j, 0 ≤ N j)
    (hcons : ∀ a b, N a = 0 → M a b * N b = 0) :
    ∀ A B, coarsenPop N f A = 0 → weightedCoarsen M N f A B * coarsenPop N f B = 0 := by
  intro A B hA
  rw [weightedCoarsen_mul_pop M N f hN A B]
  have hz : ∀ a ∈ Finset.univ.filter (fun a => f a = A), N a = 0 :=
    (Finset.sum_eq_zero_iff_of_nonneg (fun j _ => hN j)).mp (by simpa [coarsenPop] using hA)
  exact Finset.sum_eq_zero fun a ha => Finset.sum_eq_zero fun b _ => hcons a b (hz a ha)

/-- …and the converse fails: a matrix that `symmetrise` rejects can coarsen to one it
    accepts. Here an empty group reports contacts with a populated one, and merging the
    two leaves one populated group. -/
theorem weightedCoarsen_hcons_strict :
    ∃ (M : Matrix (Fin 2) (Fin 2) ℝ) (N : Fin 2 → ℝ) (f : Fin 2 → Fin 1),
      (∀ j, 0 ≤ N j) ∧ ¬ (∀ a b, N a = 0 → M a b * N b = 0) ∧
      ∀ A B, coarsenPop N f A = 0 → weightedCoarsen M N f A B * coarsenPop N f B = 0 := by
  refine ⟨fun i j => if i = 0 ∧ j = 1 then 1 else 0, fun j => if j = 0 then 0 else 1,
    fun _ => 0, ?_, ?_, ?_⟩
  · intro j; by_cases h : j = 0 <;> simp [h]
  · intro h
    have := h 0 1 (by simp)
    simp at this
  · intro A B hA
    exfalso
    have hpop : coarsenPop (fun j : Fin 2 => if j = 0 then (0 : ℝ) else 1)
        (fun _ => (0 : Fin 1)) A = 1 := by
      fin_cases A
      simp [coarsenPop, Fin.sum_univ_two]
    rw [hpop] at hA
    norm_num at hA

/-- In `MeanContacts`, converting symmetrised total contacts back is `symmetriseSafe`,
    with no hypothesis. -/
theorem ofCounts_symCounts_toCounts {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ)
    (N : Fin n → ℝ) :
    ofCounts (symCounts (toCounts M N)) N = symmetriseSafe M N := by
  funext i j
  unfold ofCounts symCounts toCounts symmetriseSafe
  split_ifs with h
  · rw [h, div_zero]
  · rw [div_div]

/-- `symmetriseSafe` is the identity on matrices that are reciprocal and in normal
    form. -/
theorem symmetriseSafe_of_reciprocal {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ)
    (N : Fin n → ℝ) (hrec : ∀ i j, M i j * N j = M j i * N i)
    (hnf : ∀ i j, N j = 0 → M i j = 0) :
    symmetriseSafe M N = M := by
  funext i j
  unfold symmetriseSafe
  split_ifs with h
  · exact (hnf i j h).symm
  · rw [← hrec i j]
    field_simp
    ring

-- ═══════════════════════════════════════════════════════════════════════════
-- `⊕`: pushforward is additive
-- ═══════════════════════════════════════════════════════════════════════════

theorem coarsenTotal_add {m n : ℕ} (A B : Matrix (Fin m) (Fin m) ℝ) (f : Fin m → Fin n) :
    coarsenTotal (A + B) f = coarsenTotal A f + coarsenTotal B f := by
  ext I J
  simp only [coarsenTotal, Matrix.add_apply, Finset.sum_add_distrib]

theorem coarsenTotal_zero {m n : ℕ} (f : Fin m → Fin n) :
    coarsenTotal (0 : Matrix (Fin m) (Fin m) ℝ) f = 0 := by
  ext I J; simp [coarsenTotal]

/-- Pushforward along `f` as an additive monoid homomorphism of total contacts. -/
noncomputable def coarsenTotalHom {m n : ℕ} (f : Fin m → Fin n) :
    Matrix (Fin m) (Fin m) ℝ →+ Matrix (Fin n) (Fin n) ℝ where
  toFun C := coarsenTotal C f
  map_zero' := coarsenTotal_zero f
  map_add' A B := coarsenTotal_add A B f

/-- Weighted coarsening is additive at a fixed population, with no hypothesis. -/
theorem weightedCoarsen_add {m n : ℕ} (A B : Matrix (Fin m) (Fin m) ℝ) (N : Fin m → ℝ)
    (f : Fin m → Fin n) :
    weightedCoarsen (A + B) N f = weightedCoarsen A N f + weightedCoarsen B N f := by
  funext I J
  simp only [weightedCoarsen, coarsenTotal, Matrix.add_apply, add_mul,
    Finset.sum_add_distrib, add_div]

-- ═══════════════════════════════════════════════════════════════════════════
-- `↑`: refinement constructs a structure with a morphism back
-- ═══════════════════════════════════════════════════════════════════════════

/-- **Refinement comes with a morphism back.** On nonnegative fine populations, `f` is a
    morphism from the structure of the refined matrix to the structure of the coarse
    matrix at the aggregated fine population. That target is the coarse matrix's own
    structure exactly when the aggregated fine population equals its population.

    `hrow` is the row half of the guard `refine` applies: an empty coarse group carries
    no contacts in its row. The column half is not needed here, because in total
    contacts a coarse entry at an empty participant group is `Mc I J · 0`; `refine`
    needs it for the matrix-level section law `weightedCoarsen_refine`. -/
noncomputable def refineHom {m n : ℕ} (Mc : Matrix (Fin n) (Fin n) ℝ) (N : Fin m → ℝ)
    (f : Fin m → Fin n) (hN : ∀ j, 0 ≤ N j)
    (hrow : ∀ I J, coarsenPop N f I = 0 → Mc I J = 0) :
    ofMean m (refineMatrix Mc N f) N ⟶ ofMean n Mc (coarsenPop N f) :=
  ⟨f, rfl, by
    funext I J
    have h := weightedCoarsen_mul_pop (refineMatrix Mc N f) N f hN I J
    -- Replace `Mc` by its normal form at the aggregated population. Refinement does not
    -- read the columns this changes, and the morphism does not see them.
    let Mc' : Matrix (Fin n) (Fin n) ℝ :=
      fun I J => if coarsenPop N f J = 0 then 0 else Mc I J
    have hrefine : refineMatrix Mc N f = refineMatrix Mc' N f := by
      funext i j
      unfold refineMatrix
      split_ifs with hj
      · rfl
      · simp only [Mc']
        rw [if_neg]
        intro hJ
        apply hj
        have hz : ∀ k ∈ Finset.univ.filter (fun k => f k = f j), N k = 0 :=
          (Finset.sum_eq_zero_iff_of_nonneg (fun k _ => hN k)).mp
            (by simpa [coarsenPop] using hJ)
        exact hz j (by simp)
    rw [hrefine, weightedCoarsen_refine Mc' N f
      (fun I J hI => by simp only [Mc']; split_ifs <;> simp [hrow I J hI])
      (fun I J hJ => by simp only [Mc']; rw [if_pos hJ])] at h
    show coarsenTotal (fun i j => refineMatrix Mc N f i j * N j) f I J =
      Mc I J * coarsenPop N f J
    rw [hrefine, ← h]
    simp only [Mc']
    split_ifs with hJ
    · rw [hJ, mul_zero, mul_zero]
    · rfl⟩

-- ═══════════════════════════════════════════════════════════════════════════
-- `⊠`: block lifts construct a structure with a morphism back
-- ═══════════════════════════════════════════════════════════════════════════

/-- `coarsenCounts` (`ConstrainedLift.lean`) is `coarsenTotal`, so the block-assembly
    results there are results about this category's pushforward. -/
theorem coarsenCounts_eq_coarsenTotal {m n : ℕ} (C : Matrix (Fin m) (Fin m) ℝ)
    (f : Fin m → Fin n) : coarsenCounts C f = coarsenTotal C f := rfl

/-- A block lift of base total contacts `C` along `f : Fin m → Fin n`: the entry for a
    fine pair `(i, j)` is the base total of its block, `C (f i) (f j)`, split by the
    kernel weight `κ (f i) (f j) i j`. `generalize` (`src/generalized_lift.jl`) is this
    construction, with `f` the projection `P × Q → P` and one kernel per base block
    over the added groups. -/
noncomputable def blockLift {m n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (κ : Fin n → Fin n → Fin m → Fin m → ℝ) (f : Fin m → Fin n) :
    Matrix (Fin m) (Fin m) ℝ :=
  fun i j => C (f i) (f j) * κ (f i) (f j) i j

/-- A block lift whose kernels sum to one over every block that carries contacts pushes
    forward to its base totals. A block with no contacts needs no condition, matching
    `generalize`'s `total == 0 && continue`, so nothing is assumed about the map. -/
theorem coarsenTotal_blockLift {m n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (κ : Fin n → Fin n → Fin m → Fin m → ℝ) (f : Fin m → Fin n)
    (hκ : ∀ I J, C I J ≠ 0 → ∑ i ∈ Finset.univ.filter (fun i => f i = I),
      ∑ j ∈ Finset.univ.filter (fun j => f j = J), κ I J i j = 1) :
    coarsenTotal (blockLift C κ f) f = C := by
  rw [← coarsenCounts_eq_coarsenTotal]
  refine coarsen_block_assembly f C (blockLift C κ f) fun I J => ?_
  calc ∑ i ∈ Finset.univ.filter (fun i => f i = I),
        ∑ j ∈ Finset.univ.filter (fun j => f j = J), blockLift C κ f i j
      = ∑ i ∈ Finset.univ.filter (fun i => f i = I),
          ∑ j ∈ Finset.univ.filter (fun j => f j = J), C I J * κ I J i j := by
        refine Finset.sum_congr rfl fun i hi => Finset.sum_congr rfl fun j hj => ?_
        rw [Finset.mem_filter] at hi hj
        show C (f i) (f j) * κ (f i) (f j) i j = _
        rw [hi.2, hj.2]
    _ = C I J * ∑ i ∈ Finset.univ.filter (fun i => f i = I),
          ∑ j ∈ Finset.univ.filter (fun j => f j = J), κ I J i j := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun i _ => (Finset.mul_sum _ _ _).symm
    _ = C I J := by
        by_cases h : C I J = 0
        · rw [h, zero_mul]
        · rw [hκ I J h, mul_one]

/-- A block lift is reciprocal (symmetric in total contacts) when its base totals are
    symmetric and each kernel is the transpose of its mirror block's. -/
theorem blockLift_symm {m n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (κ : Fin n → Fin n → Fin m → Fin m → ℝ) (f : Fin m → Fin n)
    (hC : ∀ I J, C I J = C J I) (hκ : ∀ I J i j, κ I J i j = κ J I j i) (i j : Fin m) :
    blockLift C κ f i j = blockLift C κ f j i := by
  unfold blockLift
  rw [hC (f i) (f j), hκ (f i) (f j) i j]

/-- **A normalised block lift comes with a morphism back to its base totals**: the
    projection `f` is a morphism from the lift, at a population `Nlift`, to the base
    totals at the pushforward of `Nlift`. -/
noncomputable def blockLiftHom {m n : ℕ} (C : Matrix (Fin n) (Fin n) ℝ)
    (κ : Fin n → Fin n → Fin m → Fin m → ℝ) (f : Fin m → Fin n) (Nlift : Fin m → ℝ)
    (hκ : ∀ I J, C I J ≠ 0 → ∑ i ∈ Finset.univ.filter (fun i => f i = I),
      ∑ j ∈ Finset.univ.filter (fun j => f j = J), κ I J i j = 1) :
    (⟨m, Nlift, blockLift C κ f⟩ : ContactObj) ⟶ ⟨n, coarsenPop Nlift f, C⟩ :=
  ⟨f, rfl, coarsenTotal_blockLift C κ f hκ⟩

/-- **`generalize` in `MeanContacts`** (`⊠` with a `GeneralizedLift`). Lift the
    symmetrised total contacts of a base
    matrix `M`, at the pushforward `N` of the lift population, through normalised
    kernels; convert to mean contacts at the lift population; then coarsen along `f`.
    The result is `symmetriseSafe M N`, which is `M` itself when `M` is reciprocal at `N`
    and in normal form (`symmetriseSafe_of_reciprocal`).

    `hsupp` is the condition under which Julia converts the lifted totals to mean
    contacts without throwing (`_mean_contacts_from_total_counts`): a lifted group with
    no population receives no contacts. No population sign condition is needed. -/
theorem weightedCoarsen_generalize {m n : ℕ}
    (M : Matrix (Fin n) (Fin n) ℝ) (Nlift : Fin m → ℝ)
    (κ : Fin n → Fin n → Fin m → Fin m → ℝ) (f : Fin m → Fin n)
    (hκ : ∀ I J, symCounts (toCounts M (coarsenPop Nlift f)) I J ≠ 0 →
      ∑ i ∈ Finset.univ.filter (fun i => f i = I),
      ∑ j ∈ Finset.univ.filter (fun j => f j = J), κ I J i j = 1)
    (hsupp : ∀ i j, Nlift j = 0 →
      blockLift (symCounts (toCounts M (coarsenPop Nlift f))) κ f i j = 0) :
    weightedCoarsen
        (ofCounts (blockLift (symCounts (toCounts M (coarsenPop Nlift f))) κ f) Nlift)
        Nlift f
      = symmetriseSafe M (coarsenPop Nlift f) := by
  rw [weightedCoarsen_eq_ofCounts, toCounts_ofCounts _ _ hsupp,
    coarsenTotal_blockLift _ κ f hκ, ofCounts_symCounts_toCounts]

-- ═══════════════════════════════════════════════════════════════════════════
-- Support: the zero-population guards, in total contacts
-- ═══════════════════════════════════════════════════════════════════════════

/-- Column support: nonnegative populations, and an empty participant group reports no
    contacts. These are exactly the structures of mean-contacts matrices
    (`ofMean_colSupported`, `toCounts_ofCounts`). -/
def ColSupported (A : ContactObj) : Prop :=
  (∀ j, 0 ≤ A.pop j) ∧ ∀ i j, A.pop j = 0 → A.counts i j = 0

/-- Two-sided support: nonnegative populations, and an empty group carries no contacts
    in its row or its column. On the structure of a mean-contacts matrix this is exactly
    the domain of Julia's `symmetrise` (`supported_ofMean_iff`). `refine` accepts a
    smaller set: it also requires normal form, which total contacts cannot see. -/
def Supported (A : ContactObj) : Prop :=
  (∀ j, 0 ≤ A.pop j) ∧ ∀ i j, (A.pop i = 0 ∨ A.pop j = 0) → A.counts i j = 0

private lemma fibre_zero {m n : ℕ} (N : Fin m → ℝ) (f : Fin m → Fin n)
    (hN : ∀ j, 0 ≤ N j) (J : Fin n) (hJ : coarsenPop N f J = 0) :
    ∀ j ∈ Finset.univ.filter (fun j => f j = J), N j = 0 :=
  (Finset.sum_eq_zero_iff_of_nonneg (fun j _ => hN j)).mp (by simpa [coarsenPop] using hJ)

theorem ColSupported.push {A : ContactObj} {m : ℕ} (h : ColSupported A)
    (f : Fin A.n → Fin m) : ColSupported (A.push f) := by
  refine ⟨fun J => show 0 ≤ coarsenPop A.pop f J from
    Finset.sum_nonneg (fun j _ => h.1 j), fun I J hJ => ?_⟩
  have hz := fibre_zero A.pop f h.1 J hJ
  show coarsenTotal A.counts f I J = 0
  exact Finset.sum_eq_zero fun i _ => Finset.sum_eq_zero fun j hj => h.2 i j (hz j hj)

theorem Supported.push {A : ContactObj} {m : ℕ} (h : Supported A) (f : Fin A.n → Fin m) :
    Supported (A.push f) := by
  refine ⟨fun J => show 0 ≤ coarsenPop A.pop f J from
    Finset.sum_nonneg (fun j _ => h.1 j), fun I J hIJ => ?_⟩
  show coarsenTotal A.counts f I J = 0
  rcases hIJ with hI | hJ
  · have hz := fibre_zero A.pop f h.1 I hI
    exact Finset.sum_eq_zero fun i hi =>
      Finset.sum_eq_zero fun j _ => h.2 i j (Or.inl (hz i hi))
  · have hz := fibre_zero A.pop f h.1 J hJ
    exact Finset.sum_eq_zero fun i _ =>
      Finset.sum_eq_zero fun j hj => h.2 i j (Or.inr (hz j hj))

theorem Supported.sym {A : ContactObj} (h : Supported A) : Supported A.sym := by
  refine ⟨h.1, fun i j hij => ?_⟩
  show (A.counts i j + A.counts j i) / 2 = 0
  rw [h.2 i j hij, h.2 j i hij.symm]; norm_num

/-- The structure of a mean-contacts matrix is column-supported. -/
theorem ofMean_colSupported {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ)
    (hN : ∀ j, 0 ≤ N j) : ColSupported (ofMean n M N) :=
  ⟨hN, fun i j hj => by
    change N j = 0 at hj
    show M i j * N j = 0
    rw [hj, mul_zero]⟩

/-- On the structure of a mean-contacts matrix, two-sided support is exactly
    nonnegative populations together with `hcons`, the condition under which Julia's
    `symmetrise` does not throw. -/
theorem supported_ofMean_iff {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (N : Fin n → ℝ) :
    Supported (ofMean n M N) ↔
      (∀ j, 0 ≤ N j) ∧ ∀ a b, N a = 0 → M a b * N b = 0 := by
  constructor
  · rintro ⟨hN, h⟩
    exact ⟨hN, fun a b ha => h a b (Or.inl ha)⟩
  · rintro ⟨hN, h⟩
    refine ⟨hN, fun i j hij => ?_⟩
    show M i j * N j = 0
    rcases hij with hi | hj
    · exact h i j hi
    · change N j = 0 at hj
      rw [hj, mul_zero]
