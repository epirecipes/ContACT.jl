"""
Activity refinement following the Britton-Ball social-activity construction.

The operation refines a base contact matrix on partition `P` to a product
partition `P × activity`, using observed respondent activity to split source
contact totals and an explicit mixing kernel to allocate unobserved contactee
activity.
"""

"""
    ActivityMixingKernel

Abstract supertype for blockwise activity coupling assumptions used by
`activity_refine`.
"""
abstract type ActivityMixingKernel end

"""Maximally assortative activity mixing, matching low with low, high with high, etc."""
struct AssortativeMixing <: ActivityMixingKernel end

"""Maximally disassortative activity mixing, matching low with high where possible."""
struct DisassortativeMixing <: ActivityMixingKernel end

"""Proportionate activity mixing, using the product coupling of activity marginals."""
struct ProportionateMixing <: ActivityMixingKernel end

function _activity_mixing(mixing::ActivityMixingKernel)
    mixing
end

function _activity_mixing(mixing::Symbol)
    if mixing in (:assortative, :assort, :A)
        AssortativeMixing()
    elseif mixing in (:disassortative, :disassort, :D)
        DisassortativeMixing()
    elseif mixing in (:proportionate, :proportional, :P)
        ProportionateMixing()
    else
        throw(ArgumentError("unknown activity mixing kernel :$mixing"))
    end
end

function _activity_labels(n::Integer, labels)
    if labels === nothing
        n == 2 && return ["low", "high"]
        return ["q$i" for i in 1:n]
    end
    length(labels) == n || throw(DimensionMismatch(
        "activity labels length $(length(labels)) does not match $n activity groups"))
    String.(labels)
end

"""
    ActivityRefinement(survey; n=2, cutpoints=nothing, labels=nothing,
                       mixing=:assortative, score_col=nothing)

Specification for a Britton-Ball style activity refinement.

If `cutpoints` is omitted, participants are split within each base partition
group into `n` equal-frequency activity strata by total contact count (or by
`score_col`, if supplied). Boundary participants are split fractionally, matching
the high/low convention in Britton and Ball for odd group sizes. If `cutpoints`
is supplied, participants are assigned deterministically to the intervals
`(-Inf, c₁]`, `(c₁, c₂]`, ..., `(cₖ, Inf)`.
"""
struct ActivityRefinement{K<:ActivityMixingKernel}
    survey::ContactSurvey
    n::Int
    cutpoints::Union{Nothing,Vector{Float64}}
    labels::Vector{String}
    mixing::K
    score_col::Union{Nothing,Symbol}
end

function ActivityRefinement(survey::ContactSurvey;
                            n::Integer=2,
                            cutpoints::Union{Nothing,AbstractVector{<:Real}}=nothing,
                            labels=nothing,
                            mixing::Union{Symbol,ActivityMixingKernel}=:assortative,
                            score_col::Union{Nothing,Symbol}=nothing)
    if cutpoints === nothing
        n >= 1 || throw(ArgumentError("n must be positive"))
        cps = nothing
        n_activity = Int(n)
    else
        cps = Float64.(cutpoints)
        all(isfinite, cps) || throw(ArgumentError("activity cutpoints must be finite"))
        issorted(cps) || throw(ArgumentError("activity cutpoints must be sorted"))
        allunique(cps) || throw(ArgumentError("activity cutpoints must be unique"))
        n_activity = length(cps) + 1
    end

    ActivityRefinement(
        survey,
        n_activity,
        cps,
        _activity_labels(n_activity, labels),
        _activity_mixing(mixing),
        score_col,
    )
end

"""
    activity_partition(spec::ActivityRefinement)
    activity_partition(survey, base_partition; kwargs...)

Return the activity partition for `spec`, or the refined product partition
`base_partition × activity` when called with a survey and base partition.
"""
function activity_partition(spec::ActivityRefinement)
    CategoricalPartition{:activity,Int}(collect(1:spec.n);
        participant_col=:activity, contact_col=:activity, labels=spec.labels)
end

function activity_partition(survey::ContactSurvey, base_partition::AbstractPartition; kwargs...)
    base_partition × activity_partition(ActivityRefinement(survey; kwargs...))
end

"""
    activity_mixing_plan(row_marginal, col_marginal, mixing=:assortative)

Construct a non-negative matrix with the requested row and column marginals.

Rows represent target activity strata and columns represent source activity
strata within one base contact block. Built-in kernels are `:assortative`,
`:disassortative`, and `:proportionate`.
"""
function activity_mixing_plan(row_marginal::AbstractVector{<:Real},
                              col_marginal::AbstractVector{<:Real},
                              mixing::Union{Symbol,ActivityMixingKernel}=:assortative)
    row = Float64.(row_marginal)
    col = Float64.(col_marginal)
    length(row) == length(col) || throw(DimensionMismatch(
        "row and column activity marginals must have the same length"))
    all(x -> isfinite(x) && x >= 0, row) ||
        throw(ArgumentError("row activity marginals must be finite and non-negative"))
    all(x -> isfinite(x) && x >= 0, col) ||
        throw(ArgumentError("column activity marginals must be finite and non-negative"))

    row_total = sum(row)
    col_total = sum(col)
    isapprox(row_total, col_total; rtol=1e-8, atol=1e-10) || throw(ArgumentError(
        "activity marginals must have equal totals, got $row_total and $col_total"))

    kernel = _activity_mixing(mixing)
    _activity_mixing_plan(row, col, kernel)
end

function _activity_mixing_plan(row::Vector{Float64}, col::Vector{Float64},
                               ::ProportionateMixing)
    total = (sum(row) + sum(col)) / 2
    total == 0 && return zeros(Float64, length(row), length(col))
    row * transpose(col) ./ total
end

function _activity_mixing_plan(row::Vector{Float64}, col::Vector{Float64},
                               ::AssortativeMixing)
    _ordered_activity_transport(row, col, collect(eachindex(col)))
end

function _activity_mixing_plan(row::Vector{Float64}, col::Vector{Float64},
                               ::DisassortativeMixing)
    _ordered_activity_transport(row, col, reverse(collect(eachindex(col))))
end

function _ordered_activity_transport(row::Vector{Float64}, col::Vector{Float64},
                                     col_order::Vector{Int})
    plan = zeros(Float64, length(row), length(col))
    remaining_col = copy(col)
    col_pos = 1

    for i in eachindex(row)
        remaining_row = row[i]
        while remaining_row > 1e-12 && col_pos <= length(col_order)
            j = col_order[col_pos]
            amount = min(remaining_row, remaining_col[j])
            plan[i, j] += amount
            remaining_row -= amount
            remaining_col[j] -= amount
            if remaining_col[j] <= 1e-12
                col_pos += 1
            end
        end
    end
    plan
end

"""
    activity_refine(survey, cm; kwargs...)
    activity_refine(cm, spec::ActivityRefinement)

Refine `cm` by respondent social activity following Britton and Ball.

The input matrix must use `MeanContacts` semantics and already satisfy
reciprocity in total-contact space. Apply `↔(cm)` first when starting from a raw
survey estimate. The refined matrix lives over `cm.partition × activity` and
coarsens back to `cm`.
"""
function activity_refine(survey::ContactSurvey, cm::ContactMatrix; kwargs...)
    activity_refine(cm, ActivityRefinement(survey; kwargs...))
end

function activity_refine(cm::ContactMatrix, spec::ActivityRefinement)
    cm.semantics isa MeanContacts || throw(ArgumentError(
        "activity refinement requires MeanContacts semantics"))

    base_partition = cm.partition
    n_base = n_groups(base_partition)
    n_activity = spec.n

    C = _reciprocal_total_contacts(cm)
    scores = _activity_scores(spec.survey, spec.score_col)
    base_groups = _participant_base_groups(spec.survey, base_partition)
    memberships, shares = _activity_memberships(spec, base_groups, scores, n_base)

    for g in 1:n_base
        if population(cm)[g] > 0 && sum(shares[g, :]) == 0
            throw(ArgumentError(
                "cannot refine base group $g: no survey participants assigned to this group"))
        end
    end

    alpha = _activity_source_means(spec.survey, base_partition, base_groups, memberships)
    source_marginals = _source_activity_marginals(C, shares, alpha)

    activity = activity_partition(spec)
    refined_partition = base_partition × activity
    refined_pop = _refined_activity_population(population(cm), shares)
    refined_counts = zeros(Float64, n_base * n_activity, n_base * n_activity)

    for target_base in 1:n_base
        for source_base in target_base:n_base
            row_marginal = source_marginals[source_base, target_base, :]
            col_marginal = source_marginals[target_base, source_base, :]
            plan = activity_mixing_plan(row_marginal, col_marginal, spec.mixing)

            for target_activity in 1:n_activity
                target_idx = _cartesian_index((target_base, target_activity), (n_base, n_activity))
                for source_activity in 1:n_activity
                    source_idx = _cartesian_index((source_base, source_activity), (n_base, n_activity))
                    refined_counts[target_idx, source_idx] = plan[target_activity, source_activity]
                    if target_base != source_base
                        refined_counts[source_idx, target_idx] = plan[target_activity, source_activity]
                    end
                end
            end
        end
    end

    refined_matrix = _mean_contacts_from_counts(refined_counts, refined_pop)
    ContactMatrix(refined_matrix, refined_partition, refined_pop, MeanContacts())
end

function _activity_scores(survey::ContactSurvey, score_col::Union{Nothing,Symbol})
    if score_col === nothing
        counts = Dict(id => 0.0 for id in survey.participants.part_id)
        for id in survey.contacts.part_id
            haskey(counts, id) && (counts[id] += 1.0)
        end
        return [counts[id] for id in survey.participants.part_id]
    end

    score_col in propertynames(survey.participants) ||
        throw(ArgumentError("score column :$score_col not found in participants"))
    raw_scores = survey.participants[!, score_col]
    any(ismissing, raw_scores) &&
        throw(ArgumentError("score column :$score_col contains missing values"))
    scores = Float64.(raw_scores)
    all(x -> isfinite(x) && x >= 0, scores) ||
        throw(ArgumentError("activity scores must be finite and non-negative"))
    scores
end

function _participant_base_groups(survey::ContactSurvey, base_partition::AbstractPartition)
    allunique(survey.participants.part_id) ||
        throw(ArgumentError("participant part_id values must be unique for activity refinement"))
    [assign_participant_group(base_partition, row) for row in eachrow(survey.participants)]
end

function _activity_memberships(spec::ActivityRefinement, base_groups, scores, n_base::Int)
    n_participants = length(scores)
    memberships = zeros(Float64, n_participants, spec.n)

    if spec.cutpoints === nothing
        for g in 1:n_base
            idx = findall(i -> base_groups[i] == g, eachindex(base_groups))
            isempty(idx) && continue
            sorted_idx = sort(idx, by=i -> (scores[i], i))
            m = length(sorted_idx)
            for (pos, participant_idx) in enumerate(sorted_idx)
                left = pos - 1
                right = pos
                for q in 1:spec.n
                    lower = (q - 1) * m / spec.n
                    upper = q * m / spec.n
                    overlap = min(right, upper) - max(left, lower)
                    overlap > 0 && (memberships[participant_idx, q] = overlap)
                end
            end
        end
    else
        for i in eachindex(scores)
            base_groups[i] === nothing && continue
            q = searchsortedlast(spec.cutpoints, scores[i]) + 1
            memberships[i, q] = 1.0
        end
    end

    shares = zeros(Float64, n_base, spec.n)
    for g in 1:n_base
        idx = findall(i -> base_groups[i] == g, eachindex(base_groups))
        isempty(idx) && continue
        shares[g, :] .= vec(sum(memberships[idx, :]; dims=1)) ./ length(idx)
    end
    memberships, shares
end

function _activity_source_means(survey::ContactSurvey, base_partition::AbstractPartition,
                                base_groups, memberships::AbstractMatrix{<:Real})
    n_base = n_groups(base_partition)
    n_activity = size(memberships, 2)
    denom = zeros(Float64, n_base, n_activity)
    counts = zeros(Float64, n_base, n_base, n_activity)

    for participant_idx in eachindex(base_groups)
        source_base = base_groups[participant_idx]
        source_base === nothing && continue
        for q in 1:n_activity
            denom[source_base, q] += memberships[participant_idx, q]
        end
    end

    part_id_to_idx = Dict(id => i for (i, id) in enumerate(survey.participants.part_id))
    for row in eachrow(survey.contacts)
        participant_idx = get(part_id_to_idx, row.part_id, nothing)
        participant_idx === nothing && continue
        source_base = base_groups[participant_idx]
        source_base === nothing && continue
        target_base = assign_contact_group(base_partition, row)
        target_base === nothing && continue

        for q in 1:n_activity
            counts[target_base, source_base, q] += memberships[participant_idx, q]
        end
    end

    alpha = zeros(Float64, n_base, n_base, n_activity)
    for source_base in 1:n_base
        for q in 1:n_activity
            denom[source_base, q] == 0 && continue
            alpha[:, source_base, q] .= counts[:, source_base, q] ./ denom[source_base, q]
        end
    end
    alpha
end

function _reciprocal_total_contacts(cm::ContactMatrix)
    M = matrix(cm)
    pop = population(cm)
    C = zeros(Float64, n_groups(cm), n_groups(cm))
    for j in 1:n_groups(cm)
        C[:, j] .= M[:, j] .* pop[j]
    end
    isapprox(C, transpose(C); rtol=1e-8, atol=1e-8) || throw(ArgumentError(
        "activity refinement requires reciprocal total contacts; apply ↔(cm) first"))
    C
end

function _source_activity_marginals(C::AbstractMatrix{<:Real},
                                    shares::AbstractMatrix{<:Real},
                                    alpha::Array{Float64,3})
    n_base, n_activity = size(shares)
    source_marginals = zeros(Float64, n_base, n_base, n_activity)

    for target_base in 1:n_base
        for source_base in 1:n_base
            total = C[target_base, source_base]
            total == 0 && continue

            weights = shares[source_base, :] .* alpha[target_base, source_base, :]
            denom = sum(weights)
            denom > 0 || throw(ArgumentError(
                "cannot split contacts from base group $source_base to $target_base: no observed activity-specific contacts"))

            for q in 1:n_activity
                source_marginals[target_base, source_base, q] = total * weights[q] / denom
            end
        end
    end
    source_marginals
end

function _refined_activity_population(base_population::AbstractVector{<:Real},
                                      shares::AbstractMatrix{<:Real})
    n_base, n_activity = size(shares)
    pop = zeros(Float64, n_base * n_activity)
    for g in 1:n_base
        for q in 1:n_activity
            idx = _cartesian_index((g, q), (n_base, n_activity))
            pop[idx] = base_population[g] * shares[g, q]
        end
    end
    pop
end

function _mean_contacts_from_counts(counts::AbstractMatrix{<:Real},
                                    pop::AbstractVector{<:Real})
    M = zeros(Float64, size(counts))
    for j in eachindex(pop)
        if pop[j] == 0
            all(iszero, counts[:, j]) || throw(ArgumentError(
                "cannot convert refined counts for zero-population activity group $j"))
        else
            M[:, j] .= counts[:, j] ./ pop[j]
        end
    end
    M
end
