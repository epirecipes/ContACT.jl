"""
Matrix construction: ContactSurvey → ContactMatrix.

At a fixed partition a survey is summarised by two sufficient statistics: the
weighted contact counts `C[i,j]`, from participant group `j` into contactee group
`i`, and the weighted participant mass `W_j` per group. The matrix is the
normalisation `M[i,j] = C[i,j]/W_j`.

Both statistics push forward along a partition map, but the matrix is normalised by
participant mass while `coarsen` weights by population, so the two only agree under
the calibration condition documented on `compute_matrix`.
"""

"""
    compute_matrix(survey::ContactSurvey, partition::AbstractPartition;
                   population=nothing, weights=nothing)

Compute a contact matrix from survey data over a fixed finite partition.

Entry `M[i,j]` is the weighted number of contacts reported into group `i` by
participants in group `j`, divided by the weighted participant mass `W_j` of group
`j`. A group with no participants gets a zero column.

# Relation to coarsening

Building at a fine partition and coarsening along `f` is guaranteed to give the same
*matrix* as building directly at the coarse partition whenever, within each coarse
participant/source fibre `J`, the attached population is proportional to participant
mass:

    N_j = λ_J · W_j        for all j ∈ f⁻¹(J)

with `λ_J > 0`, required only on fibres carrying participant mass — on the others,
`N` is unconstrained. The constant cancels in the derivation, but it must be
nonzero; on ContACT's nonnegative population domain that means strictly positive.
The condition is per-fibre, so proportionality must hold separately in each one.

This condition is sufficient for every survey and is required for agreement for
arbitrary contact-count data. A particular survey may agree without it — for example,
when all fine participant groups in a fibre have the same coarse contact profile —
but the operations do not guarantee that. A fibre with no participant mass
contributes zero source columns on both routes, so the matrix values agree there for
any `N`.

Matrix agreement does not by itself imply agreement of the attached populations. To
guarantee agreement of the full `ContactMatrix` objects under an external `N`, the
direct coarse build should also be given the pushed-forward population `f_*N`.

With `population=nothing`, the condition holds automatically because the weighted
participant mass `W` is itself attached as the population. With an external
`population`, it becomes a calibration (post-stratification) condition on the survey
weights. Supplying `weights` alone does not establish it: when weighted sample
composition differs from the target population within a fibre, the two routes
generally estimate different quantities and need not agree.

# Arguments
- `survey`: a `ContactSurvey` with participant and contact data
- `partition`: the finite partition to use

# Keyword Arguments
- `population`: population per partition group to attach to the resulting matrix.
  If `nothing`, uses the weighted participant mass computed from the survey.
- `weights`: participant column containing sampling weights, or `nothing` for equal
  participant weights

# Returns
A `ContactMatrix{Float64, MeanContacts}` giving the mean number of contacts reported
by participants in each source group with contacts in each target group.

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
    dropped_missing_contact_group = 0
    dropped_unknown_participant = 0
    dropped_missing_participant_group = 0

    for (k, cnt_grp) in enumerate(cnt_groups)
        if cnt_grp === nothing
            dropped_missing_contact_group += 1
            continue
        end
        part_id = survey.contacts.part_id[k]
        idx = get(part_id_to_idx, part_id, nothing)
        if idx === nothing
            dropped_unknown_participant += 1
            continue
        end
        part_grp = part_groups[idx]
        if part_grp === nothing
            dropped_missing_participant_group += 1
            continue
        end
        contact_counts[cnt_grp, part_grp] += w[idx]
    end
    _warn_dropped_contacts("compute_matrix";
        missing_contact_group=dropped_missing_contact_group,
        unknown_participant=dropped_unknown_participant,
        missing_participant_group=dropped_missing_participant_group)

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

function _warn_dropped_contacts(context::AbstractString;
                                missing_contact_group::Int=0,
                                unknown_participant::Int=0,
                                missing_participant_group::Int=0)
    dropped = missing_contact_group + unknown_participant + missing_participant_group
    dropped == 0 && return nothing

    reasons = String[]
    missing_contact_group == 0 || push!(reasons,
        "$missing_contact_group with missing/unmapped contact group")
    unknown_participant == 0 || push!(reasons,
        "$unknown_participant with unknown participant id")
    missing_participant_group == 0 || push!(reasons,
        "$missing_participant_group with missing/unmapped participant group")
    detail = join(reasons, "; ")
    @warn "$context dropped $dropped contact(s): $detail"
    nothing
end
