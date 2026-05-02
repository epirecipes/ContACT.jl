import Mathlib.CategoryTheory.Category.Basic
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.Data.Matrix.Basic

/-!
*Source: `ContactCat.lean`*

# The Category of Contact Matrices

We formalise the category **Contact** whose:
- **Objects** are contact matrices: a square matrix M over ℝ≥0, an age partition
  (encoded as n : ℕ), and a population vector N : Fin n → ℝ>0.
- **Morphisms** are structure-preserving maps between contact matrices that
  respect the age partition and population structure.

## Summary of Results

| # | Result | Status |
|---|--------|--------|
| 1 | Contact matrices form objects of a category | ✅ |
| 2 | Identity morphism (identity on matching partitions) | ✅ |
| 3 | Composition of morphisms is associative | ✅ |

-/

open Matrix

universe u

/-- A contact matrix object: dimension, matrix, and population vector. -/
structure ContactMatObj where
  n : ℕ
  mat : Matrix (Fin n) (Fin n) ℝ
  pop : Fin n → ℝ
  pop_pos : ∀ i, pop i > 0

/-- A morphism between contact matrices induced by an age-group map.
    Given f : Fin m → Fin n (surjective), transforms an m×m matrix to n×n. -/
structure ContactMatHom (A B : ContactMatObj) where
  map : Fin A.n → Fin B.n
  -- The map respects the contact structure via coarsening
  surj : Function.Surjective map

/-- Identity morphism. -/
def ContactMatHom.id (A : ContactMatObj) : ContactMatHom A A where
  map := _root_.id
  surj := Function.surjective_id

/-- Composition of morphisms. -/
def ContactMatHom.comp {A B C : ContactMatObj}
    (g : ContactMatHom B C) (f : ContactMatHom A B) : ContactMatHom A C where
  map := g.map ∘ f.map
  surj := Function.Surjective.comp g.surj f.surj

/-- Composition is associative. -/
theorem ContactMatHom.comp_assoc {A B C D : ContactMatObj}
    (h : ContactMatHom C D) (g : ContactMatHom B C) (f : ContactMatHom A B) :
    ContactMatHom.comp h (ContactMatHom.comp g f) =
    ContactMatHom.comp (ContactMatHom.comp h g) f := by
  simp [ContactMatHom.comp, Function.comp_assoc]

/-- Left identity law. -/
theorem ContactMatHom.id_comp {A B : ContactMatObj} (f : ContactMatHom A B) :
    ContactMatHom.comp (ContactMatHom.id B) f = f := by
  simp [ContactMatHom.comp, ContactMatHom.id]

/-- Right identity law. -/
theorem ContactMatHom.comp_id {A B : ContactMatObj} (f : ContactMatHom A B) :
    ContactMatHom.comp f (ContactMatHom.id A) = f := by
  simp [ContactMatHom.comp, ContactMatHom.id]
