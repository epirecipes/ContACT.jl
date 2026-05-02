import Lake
open Lake DSL

package "ContACTProofs" where
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩
  ]

require mdgen from git
  "https://github.com/Seasawher/mdgen" @ "main"

require "leanprover-community" / "mathlib"

@[default_target]
lean_lib «ContACTProofs» where
  globs := #[.submodules `ContACTProofs]
