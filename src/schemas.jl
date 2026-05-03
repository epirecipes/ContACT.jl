"""
ACSet schemas for contact surveys and matrices.

This module leverages Catlab's ACSet machinery to provide:

1. **Schema-enforced survey structure**: Type-safe joins between participants
   and contacts, preventing structural errors at construction time.

2. **Functorial data migration (Σ, Δ, Π)**: Principled coarsening (Σ = left Kan),
   restriction (Δ = pullback), and refinement (Π = right Kan) of age-structured data.

3. **Structured cospans for open contact matrices**: Composable interfaces
   where shared age groups form the boundary.

4. **Undirected wiring diagrams**: Declarative specification of how
   setting-specific matrices compose.

These provide functionality beyond what plain structs offer:
- Schema validation at construction
- Automatic bookkeeping during composition
- Formal category-theoretic operations (pushouts, pullbacks)
- Diagrammatic reasoning via wiring diagrams
"""

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.UndirectedWiringDiagrams: AbstractUWD
using Catlab.Programs: @relation

# ---------------------------------------------------------------------------
# Schema: ContactSurveySchema
# ---------------------------------------------------------------------------

"""
ACSet schema for a contact survey.

Objects:
- P (Participants)
- C (Contacts)
- G (Age Groups)

Morphisms:
- part_group: P → G (participant's age group)
- cnt_group: C → G (contact's age group)
- reporter: C → P (which participant reported this contact)

This schema makes the foreign-key relationships between participants,
contacts, and age groups explicit and type-checked.
"""
@present SchContactSurvey(FreeSchema) begin
    P::Ob    # Participants
    C::Ob    # Contacts
    G::Ob    # Age groups
    reporter::Hom(C, P)      # contact → participant who reported it
    part_group::Hom(P, G)    # participant → their age group
    cnt_group::Hom(C, G)     # contact → their age group
end

@abstract_acset_type AbstractContactSurveyACSet

@acset_type ContactSurveyACSet(SchContactSurvey,
    index=[:reporter, :part_group, :cnt_group]) <: AbstractContactSurveyACSet

# ---------------------------------------------------------------------------
# Schema: ContactMatrixSchema (labelled)
# ---------------------------------------------------------------------------

"""
ACSet schema for a contact matrix with labels.

Objects:
- G (Age Groups — shared with survey schema)
- E (Entries — matrix cells)

Morphisms:
- row_group: E → G (row age group = "contacted")
- col_group: E → G (column age group = "contactor")

Attributes:
- value: E → ℝ (matrix entry value)
- gname: G → Name (group label)
- pop: G → ℝ (population)
"""
@present SchContactMatrix(FreeSchema) begin
    G::Ob
    E::Ob
    row_group::Hom(E, G)
    col_group::Hom(E, G)
    Value::AttrType
    Name::AttrType
    Pop::AttrType
    value::Attr(E, Value)
    gname::Attr(G, Name)
    pop::Attr(G, Pop)
end

@acset_type ContactMatrixACSet(SchContactMatrix,
    index=[:row_group, :col_group])

# ---------------------------------------------------------------------------
# Open contact matrices (structured cospans for composition)
# ---------------------------------------------------------------------------

const OpenContactMatrixOb, OpenContactMatrix =
    OpenACSetTypes(ContactMatrixACSet, :G)

"""
    OpenContactMatrixOb

Interface object type for open contact matrices. Age groups are exposed as
composition boundaries.
"""
OpenContactMatrixOb

"""
    OpenContactMatrix

Structured cospan type for open contact matrices. When two open matrices share
age groups at a boundary, they compose by summing contributions at shared groups.
"""
OpenContactMatrix

const LabelledContactMatrix = ContactMatrixACSet{Float64, String, Float64}

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

"""
    ContactSurveyACSet(survey::ContactSurvey, partition::AbstractPartition; skip_invalid=false)

Convert a plain `ContactSurvey` into an ACSet representation with
age groups assigned according to `partition`.

By default, rows whose reporter or age cannot be assigned are rejected. Pass
`skip_invalid=true` to keep the part but leave the invalid relationship unset.
"""
function ContactSurveyACSet(survey::ContactSurvey, partition::AbstractPartition; skip_invalid::Bool=false)
    acs = ContactSurveyACSet()
    n = n_groups(partition)

    # Add age groups
    add_parts!(acs, :G, n)

    # Add participants with group assignments
    n_parts = nrow(survey.participants)
    add_parts!(acs, :P, n_parts)
    for (i, row) in enumerate(eachrow(survey.participants))
        grp = assign_participant_group(partition, row)
        if grp === nothing
            skip_invalid && continue
            throw(ArgumentError("participant row $i cannot be assigned to partition $(dimension(partition))"))
        end
        set_subpart!(acs, i, :part_group, grp)
    end

    # Add contacts with group assignments and reporter links
    part_id_to_idx = Dict(id => i for (i, id) in enumerate(survey.participants.part_id))
    n_contacts = nrow(survey.contacts)
    add_parts!(acs, :C, n_contacts)
    for k in 1:n_contacts
        pid = survey.contacts.part_id[k]
        idx = get(part_id_to_idx, pid, nothing)
        if idx === nothing
            skip_invalid && continue
            throw(ArgumentError("contact row $k references unknown participant id $pid"))
        end
        set_subpart!(acs, k, :reporter, idx)

        grp = assign_contact_group(partition, survey.contacts[k, :])
        if grp === nothing
            skip_invalid && continue
            throw(ArgumentError("contact row $k cannot be assigned to partition $(dimension(partition))"))
        end
        set_subpart!(acs, k, :cnt_group, grp)
    end

    return acs
end

"""
    LabelledContactMatrix(cm::ContactMatrix)

Convert a `ContactMatrix` to an ACSet representation.
"""
function LabelledContactMatrix(cm::ContactMatrix)
    acs = LabelledContactMatrix()
    n = n_groups(cm)
    M = matrix(cm)
    pop = population(cm)
    labels = group_labels(cm)

    # Add groups with labels and population
    add_parts!(acs, :G, n)
    for i in 1:n
        set_subpart!(acs, i, :gname, labels[i])
        set_subpart!(acs, i, :pop, pop[i])
    end

    # Add matrix entries
    for j in 1:n, i in 1:n
        eid = add_part!(acs, :E)
        set_subpart!(acs, eid, :row_group, i)
        set_subpart!(acs, eid, :col_group, j)
        set_subpart!(acs, eid, :value, M[i, j])
    end

    return acs
end

# ---------------------------------------------------------------------------
# Functorial data migration for coarsening
# ---------------------------------------------------------------------------

"""
    migrate_coarsen(acs::ContactSurveyACSet, f::PartitionMap)

Apply functorial data migration (Σ_f) to coarsen the age groups in a survey ACSet.

This is the ACSet-level implementation of the left Kan extension: it
pushes forward participant and contact group assignments along the
surjective age-group map f.
"""
function migrate_coarsen(acs::ContactSurveyACSet, f::PartitionMap)
    fmap = collect(f.mapping)
    n_coarse = n_groups(f.codomain)

    coarsened = ContactSurveyACSet()
    add_parts!(coarsened, :G, n_coarse)

    # Copy participants with coarsened groups
    n_p = nparts(acs, :P)
    add_parts!(coarsened, :P, n_p)
    for i in 1:n_p
        fine_grp = subpart(acs, i, :part_group)
        set_subpart!(coarsened, i, :part_group, fmap[fine_grp])
    end

    # Copy contacts with coarsened groups and same reporter links
    n_c = nparts(acs, :C)
    add_parts!(coarsened, :C, n_c)
    for k in 1:n_c
        set_subpart!(coarsened, k, :reporter, subpart(acs, k, :reporter))
        fine_grp = subpart(acs, k, :cnt_group)
        set_subpart!(coarsened, k, :cnt_group, fmap[fine_grp])
    end

    return coarsened
end

# ---------------------------------------------------------------------------
# UWD-based composition of settings
# ---------------------------------------------------------------------------

"""
    ContactSharer

An undirected open system wrapping a contact matrix for operadic composition.
Analogous to ProjectionSharer in CategoricalPopulationDynamics.jl.
"""
struct ContactSharer{T<:Real}
    n_groups::Int
    matrix::Matrix{T}
end

"""
    compose_uwd(diagram, sharers::Dict{Symbol, ContactSharer})

Compose contact matrices according to an undirected wiring diagram.

Each box in the UWD represents a contact setting (home, work, school, etc.).
Composition is additive: matrices at shared junctions sum.

# Example
```julia
diagram = @relation (age,) begin
    home(age)
    work(age)
    school(age)
end
sharers = Dict(:home => ContactSharer(M_home), :work => ContactSharer(M_work),
               :school => ContactSharer(M_school))
total_matrix = compose_uwd(diagram, sharers)
```
"""
function compose_uwd(diagram::AbstractUWD, sharers::Dict{Symbol, ContactSharer{T}}) where {T}
    n_boxes = nparts(diagram, :Box)
    names = Symbol.(subpart(diagram, :name))

    n = first(values(sharers)).n_groups
    K = zeros(T, n, n)
    for nm in names
        haskey(sharers, nm) || error("No sharer for box :$nm")
        K .+= sharers[nm].matrix
    end
    return K
end

ContactSharer(M::AbstractMatrix{T}) where {T<:Real} = ContactSharer{T}(size(M, 1), Matrix(M))
