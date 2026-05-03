# Composition & Stratification

## Additive Composition (⊕)

Contact matrices from different settings (home, work, school) form a
**commutative monoid** under addition:

```@example composition
using ContACT
using LinearAlgebra

partition = AgePartition([0, 5, 18, 65])
uk_pop = [6987.0, 11537.0, 35854.0, 9492.0]

M_home = [1.8 0.3 0.4 0.1;
          0.3 0.6 0.2 0.1;
          0.4 0.2 1.0 0.3;
          0.1 0.1 0.3 0.8]

M_work = [0.0 0.0 0.0 0.0;
          0.0 0.1 1.5 0.2;
          0.0 1.5 2.5 0.4;
          0.0 0.2 0.4 0.3]

M_school = [3.5 2.0 0.1 0.0;
            2.0 5.0 0.3 0.0;
            0.1 0.3 0.2 0.0;
            0.0 0.0 0.0 0.0]

M_other = [0.5 0.3 0.4 0.2;
           0.3 0.5 0.8 0.3;
           0.4 0.8 1.2 0.5;
           0.2 0.3 0.5 0.4]

cm_home   = ContactMatrix(M_home, partition, uk_pop)
cm_work   = ContactMatrix(M_work, partition, uk_pop)
cm_school = ContactMatrix(M_school, partition, uk_pop)
cm_other  = ContactMatrix(M_other, partition, uk_pop)

# Compose with ⊕
cm_total = cm_home ⊕ cm_work ⊕ cm_school ⊕ cm_other
round.(matrix(cm_total); digits=2)
```

### Monoid Laws

```@example composition
# Associativity
lhs = (cm_home ⊕ cm_work) ⊕ cm_school
rhs = cm_home ⊕ (cm_work ⊕ cm_school)
println("Associative: $(matrix(lhs) ≈ matrix(rhs))")

# Commutativity
println("Commutative: $(matrix(cm_home ⊕ cm_work) ≈ matrix(cm_work ⊕ cm_home))")

# Identity
zero_cm = ContactMatrix(zeros(4, 4), partition, uk_pop)
println("Identity: $(matrix(cm_home ⊕ zero_cm) ≈ matrix(cm_home))")
nothing # hide
```

## Intervention Modelling

The compositional structure enables scenario analysis by modifying individual
settings:

```@example composition
# School closure
cm_closure = cm_home ⊕ cm_work ⊕ cm_other

# Work-from-home (70% reduction)
cm_wfh = ContactMatrix(0.3 .* M_work, partition, uk_pop)
cm_lockdown = cm_home ⊕ cm_wfh ⊕ cm_other

println("Spectral radius ρ(M) ∝ R₀:")
println("  Normal:         $(round(ρ(cm_total); digits=2))")
println("  School closure: $(round(ρ(cm_closure); digits=2))")
println("  Lockdown:       $(round(ρ(cm_lockdown); digits=2))")
println("  Home only:      $(round(ρ(cm_home); digits=2))")
nothing # hide
```

## Spatial Stratification (⊗)

The Kronecker product with a coupling matrix creates multi-region contact
patterns:

```@example composition
# Three regions with mostly-local mixing
coupling = [0.8 0.15 0.05;
            0.15 0.7 0.15;
            0.05 0.15 0.8]

cm_spatial = cm_total ⊗ coupling
println("Stratified: $(n_groups(cm_spatial)) groups (3 regions × 4 ages)")
println("ρ(M_spatial): $(round(ρ(cm_spatial); digits=2))")
nothing # hide
```

The block structure of the stratified matrix:

```@example composition
M_spatial = matrix(cm_spatial)
println("Diagonal block (region 1 local):")
display(round.(M_spatial[1:4, 1:4]; digits=2))
println("\nOff-diagonal block (region 1 → region 2):")
display(round.(M_spatial[1:4, 5:8]; digits=2))
nothing # hide
```

## Composing Operations

Operations chain naturally:

```@example composition
# Coarsen then stratify
cm_coarse = cm_total ↓ AgePartition([0, 18])
cm_2reg = cm_coarse ⊗ [0.9 0.1; 0.1 0.9]
println("Pipeline: 4-group → 2-group → 2-region = $(n_groups(cm_2reg)) groups")
round.(matrix(cm_2reg); digits=2)
```
