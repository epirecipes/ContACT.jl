"""
The functor: ContactSurvey → ContactMatrix.

This implements the restricted functor from the subcategory of surveys
(with fixed age partition) to the category of contact matrices.

Categorically: compute_matrix is a functor F: Survey(P) → Contact
where Survey(P) is the subcategory of surveys equipped with a fixed
age partition P.
"""

"""
    compute_matrix(survey::ContactSurvey, partition::AgePartition;
                   population=nothing, weights=nothing)

Compute a contact matrix from survey data with a fixed age partition.

This is the central functor of ContACT.jl: it maps individual-level
contact data to an age-structured contact matrix.

# Arguments
- `survey`: a `ContactSurvey` with participant and contact data
- `partition`: the age discretisation to use

# Keyword Arguments
- `population`: population vector per age group (for weighting/symmetrisation).
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
function compute_matrix(survey::ContactSurvey, partition::AgePartition;
                        population::Union{Nothing, AbstractVector{<:Real}}=nothing,
                        weights::Union{Nothing, Symbol}=nothing)
    limits = age_limits(partition)
    n = n_groups(partition)

    # Assign age groups to participants
    part_groups = _assign_age_group.(survey.participants.part_age, Ref(limits))
    # Assign age groups to contacts
    cnt_groups = _assign_age_group.(survey.contacts.cnt_age, Ref(limits))

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

    # Count weighted contacts per (participant_group, contact_group) pair
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

# Assign a single age to an age group index (1-based), or nothing if missing.
function _assign_age_group(age, limits::Vector{Float64})
    age === missing && return nothing
    a = Float64(age)
    isfinite(a) || return nothing
    # Find the last limit ≤ age
    grp = searchsortedlast(limits, a)
    grp == 0 && return nothing  # age below minimum limit
    return grp
end
