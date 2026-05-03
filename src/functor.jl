"""
The functor: ContactSurvey → ContactMatrix.

This implements the restricted functor from the subcategory of surveys
(with fixed partition) to the category of contact matrices.

Categorically: compute_matrix is a functor F: Survey(P) → Contact
where Survey(P) is the subcategory of surveys equipped with a fixed
finite partition P.
"""

"""
    compute_matrix(survey::ContactSurvey, partition::AbstractPartition;
                   population=nothing, weights=nothing)

Compute a contact matrix from survey data with a fixed finite partition.

This is the central functor of ContACT.jl: it maps individual-level
contact data to a structured contact matrix.

# Arguments
- `survey`: a `ContactSurvey` with participant and contact data
- `partition`: the finite partition to use

# Keyword Arguments
- `population`: population vector per group (for weighting/symmetrisation).
  If `nothing`, uses the empirical participant distribution.
- `weights`: column name in participants to use as sampling weights, or `nothing`

# Returns
A `ContactMatrix{Float64, MeanContacts}` giving the mean number of contacts
reported by participants in each age group with contacts in each other group.

# Example
```julia
partition = AgePartition([0, 5, 18, 65])
cm = compute_matrix(survey, partition)
```
"""
function compute_matrix(survey::ContactSurvey, partition::AbstractPartition;
                        population::Union{Nothing, AbstractVector{<:Real}}=nothing,
                        weights::Union{Nothing, Symbol}=nothing)
    n = n_groups(partition)

    part_groups = [assign_participant_group(partition, row) for row in eachrow(survey.participants)]
    cnt_groups = [assign_contact_group(partition, row) for row in eachrow(survey.contacts)]

    # Build weight vector
    if weights === nothing
        w = ones(Float64, nrow(survey.participants))
    else
        weights in propertynames(survey.participants) ||
            throw(ArgumentError("weight column :$weights not found in participants"))
        raw_weights = survey.participants[!, weights]
        any(ismissing, raw_weights) &&
            throw(ArgumentError("weight column :$weights contains missing values"))
        w = Float64.(raw_weights)
        all(x -> isfinite(x) && x >= 0, w) ||
            throw(ArgumentError("weight column :$weights must be finite and non-negative"))
    end

    # Count weighted contacts per (participant group, contact group) pair
    contact_counts = zeros(Float64, n, n)
    participant_weights = zeros(Float64, n)

    # Map part_id to row index for O(1) lookup
    part_id_to_idx = Dict(id => i for (i, id) in enumerate(survey.participants.part_id))

    for (k, cnt_grp) in enumerate(cnt_groups)
        cnt_grp === nothing && continue  # skip contacts with missing age
        part_id = survey.contacts.part_id[k]
        idx = get(part_id_to_idx, part_id, nothing)
        idx === nothing && continue
        part_grp = part_groups[idx]
        part_grp === nothing && continue
        contact_counts[cnt_grp, part_grp] += w[idx]
    end

    # Sum weights per participant group
    for (i, grp) in enumerate(part_groups)
        grp === nothing && continue
        participant_weights[grp] += w[i]
    end

    # Normalise: mean contacts = total contacts / number of participants
    M = zeros(Float64, n, n)
    for j in 1:n
        if participant_weights[j] > 0
            M[:, j] .= contact_counts[:, j] ./ participant_weights[j]
        end
    end

    # Population: use empirical or provided
    pop = if population !== nothing
        Float64.(population)
    else
        participant_weights
    end

    ContactMatrix(M, partition, pop, MeanContacts())
end
