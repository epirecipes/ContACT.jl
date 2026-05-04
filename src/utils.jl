"""
Utility functions for ContACT.jl.
"""

"""
    to_per_capita(cm::ContactMatrix)

Convert a contact matrix from mean contacts to per-capita rates.
Divides each column by the population of the contactor age group.
"""
function to_per_capita(cm::ContactMatrix)
    M = matrix(cm)
    pop = population(cm)
    n = n_groups(cm)

    M_pc = zeros(Float64, n, n)
    for j in 1:n
        if pop[j] == 0
            all(iszero, M[:, j]) || throw(ArgumentError(
                "cannot convert column $j to per-capita rates: zero population with nonzero contacts"))
        else
            for i in 1:n
                M_pc[i, j] = M[i, j] / pop[j]
            end
        end
    end
    ContactMatrix(M_pc, cm.partition, pop, PerCapitaRate())
end

"""
    to_counts(cm::ContactMatrix)

Convert a contact matrix from mean contacts to total contact counts.
Multiplies each column by the population of the participant age group.
"""
function to_counts(cm::ContactMatrix)
    M = matrix(cm)
    pop = population(cm)
    n = n_groups(cm)

    M_counts = zeros(Float64, n, n)
    for j in 1:n
        M_counts[:, j] .= M[:, j] .* pop[j]
    end
    ContactMatrix(M_counts, cm.partition, pop, ContactCounts())
end

"""
    spectral_radius(cm::ContactMatrix)

Compute the dominant eigenvalue (spectral radius) of the contact matrix.
"""
function spectral_radius(cm::ContactMatrix)
    maximum(abs.(eigvals(matrix(cm))))
end

"""
    next_generation_matrix(cm; transmissibility=1, recovery_rate=1)

Return the next-generation matrix associated with a contact matrix under
frequency-dependent transmission.

For a `MeanContacts` matrix `M`, entry `(i, j)` is
`transmissibility / recovery_rate * M[i, j] * N[i] / N[j]`, the expected new
infections in target group `i` caused by one infectious individual in source
group `j` in an otherwise susceptible population.
"""
function next_generation_matrix(cm::ContactMatrix;
                                transmissibility::Real=1,
                                recovery_rate::Real=1)
    cm.semantics isa MeanContacts || throw(ArgumentError(
        "next_generation_matrix requires MeanContacts semantics"))
    isfinite(transmissibility) && transmissibility >= 0 ||
        throw(ArgumentError("transmissibility must be finite and non-negative"))
    isfinite(recovery_rate) && recovery_rate > 0 ||
        throw(ArgumentError("recovery_rate must be finite and positive"))

    M = matrix(cm)
    pop = population(cm)
    n = n_groups(cm)
    K = zeros(Float64, n, n)
    scale = Float64(transmissibility) / Float64(recovery_rate)
    for source in 1:n
        if pop[source] == 0
            all(iszero, M[:, source]) || throw(ArgumentError(
                "cannot build next-generation matrix: source group $source has zero population with nonzero contacts"))
            continue
        end
        for target in 1:n
            K[target, source] = scale * M[target, source] * pop[target] / pop[source]
        end
    end
    K
end

"""
    basic_reproduction_number(cm; transmissibility=1, recovery_rate=1)
    R0(cm; transmissibility=1, recovery_rate=1)
    R₀(cm; transmissibility=1, recovery_rate=1)

Compute the spectral radius of the next-generation matrix.
"""
function basic_reproduction_number(cm::ContactMatrix; kwargs...)
    maximum(abs.(eigvals(next_generation_matrix(cm; kwargs...))))
end

"""
    R0(cm; transmissibility=1, recovery_rate=1)

ASCII alias for `basic_reproduction_number`.
"""
R0(cm::ContactMatrix; kwargs...) = basic_reproduction_number(cm; kwargs...)

"""
    R₀(cm; transmissibility=1, recovery_rate=1)

Unicode alias for `basic_reproduction_number`.
"""
R₀(cm::ContactMatrix; kwargs...) = basic_reproduction_number(cm; kwargs...)

"""
    calibrate_transmissibility(cm, target_R0; recovery_rate=1)

Return the transmissibility parameter that gives `target_R0` for `cm` under the
next-generation matrix convention used by `R0`.
"""
function calibrate_transmissibility(cm::ContactMatrix, target_R0::Real;
                                    recovery_rate::Real=1)
    isfinite(target_R0) && target_R0 >= 0 ||
        throw(ArgumentError("target_R0 must be finite and non-negative"))
    isfinite(recovery_rate) && recovery_rate > 0 ||
        throw(ArgumentError("recovery_rate must be finite and positive"))

    baseline = basic_reproduction_number(cm; transmissibility=1, recovery_rate=1)
    baseline > 0 || throw(ArgumentError(
        "cannot calibrate transmissibility for a matrix with zero epidemic threshold"))
    Float64(target_R0) * Float64(recovery_rate) / baseline
end

"""
    marginal_matrix(cm, partition_or_dimension)

Return the contact matrix marginalized to a product factor or subset of factors.

`partition_or_dimension` can be an explicit `AbstractPartition`, a dimension
symbol such as `:sep`, or a tuple of symbols such as `(:sep, :edu)`.
"""
marginal_matrix(cm::ContactMatrix, partition::AbstractPartition) = coarsen(cm, partition)

function marginal_matrix(cm::ContactMatrix, dim::Symbol)
    marginal_matrix(cm, _partition_for_dimension(cm.partition, dim))
end

function marginal_matrix(cm::ContactMatrix, dims::Tuple)
    all(d -> d isa Symbol, dims) ||
        throw(ArgumentError("marginal dimensions must be symbols"))
    factors = Tuple(_partition_for_dimension(cm.partition, d) for d in dims)
    marginal_matrix(cm, ProductPartition(factors))
end

"""
    assortativity_index(cm, partition_or_dimension)
    assortativity_index(cm)

Compute the trace of the row-normalized marginal matrix.

For a two-group SEP or education marginal this matches the CoMix/SEP summary:
the fraction of contacts made within each group, summed across groups.
"""
assortativity_index(cm::ContactMatrix, partition_or_dimension) =
    assortativity_index(marginal_matrix(cm, partition_or_dimension))

function assortativity_index(cm::ContactMatrix)
    M = matrix(cm)
    total = 0.0
    for i in 1:n_groups(cm)
        row_total = sum(M[i, :])
        row_total > 0 || throw(ArgumentError(
            "cannot compute assortativity index: row $i has zero total contacts"))
        total += M[i, i] / row_total
    end
    total
end

"""
    type_reproduction_number(cm, target_groups; transmissibility=1, recovery_rate=1)

Compute the type-reproduction number for a targeted set of groups using the
next-generation matrix.

The target groups can be group indices, a Boolean mask, or group labels. The
returned value is the spectral radius of the target-to-target next-generation
operator after accounting for infections that pass through non-target groups.
"""
function type_reproduction_number(cm::ContactMatrix, target_groups;
                                  transmissibility::Real=1,
                                  recovery_rate::Real=1)
    K = next_generation_matrix(cm;
        transmissibility=transmissibility,
        recovery_rate=recovery_rate,
    )
    target = _group_indices(cm, target_groups)
    other = setdiff(collect(1:n_groups(cm)), target)
    isempty(target) && throw(ArgumentError("target_groups must not be empty"))

    if isempty(other)
        return maximum(abs.(eigvals(K[target, target])))
    end

    Ktt = K[target, target]
    Kto = K[target, other]
    Kot = K[other, target]
    Koo = K[other, other]
    bridge = I - Koo
    H = Ktt + Kto * (bridge \ Kot)
    maximum(abs.(eigvals(H)))
end

"""
    control_threshold(Tg)

Return the targeted-control threshold `max(0, 1 - 1/Tg)` for a type
reproduction number.
"""
function control_threshold(Tg::Real)
    isfinite(Tg) && Tg >= 0 ||
        throw(ArgumentError("type reproduction number must be finite and non-negative"))
    Tg == 0 && return 0.0
    max(0.0, 1.0 - inv(Float64(Tg)))
end

"""
    control_effort(cm, target_groups, threshold)
    control_effort(cm, target_groups; transmissibility=1, recovery_rate=1)

Compute the population-weighted fraction of the whole population controlled by a
targeted strategy. If `threshold` is omitted, it is derived from
`type_reproduction_number`.
"""
function control_effort(cm::ContactMatrix, target_groups, threshold::Real)
    isfinite(threshold) && threshold >= 0 ||
        throw(ArgumentError("control threshold must be finite and non-negative"))
    total_population = sum(population(cm))
    total_population > 0 ||
        throw(ArgumentError("cannot compute control effort for zero total population"))
    target = _group_indices(cm, target_groups)
    threshold * sum(population(cm)[target]) / total_population
end

function control_effort(cm::ContactMatrix, target_groups;
                        transmissibility::Real=1,
                        recovery_rate::Real=1)
    Tg = type_reproduction_number(cm, target_groups;
        transmissibility=transmissibility,
        recovery_rate=recovery_rate,
    )
    control_effort(cm, target_groups, control_threshold(Tg))
end

function _partition_for_dimension(partition::AbstractPartition, dim::Symbol)
    dimension(partition) == dim || throw(ArgumentError(
        "partition dimension $(dimension(partition)) does not match :$dim"))
    partition
end

function _partition_for_dimension(partition::ProductPartition, dim::Symbol)
    matches = [factor for factor in partition.factors if dimension(factor) == dim]
    length(matches) == 1 || throw(ArgumentError(
        "dimension :$dim must match exactly one product factor"))
    only(matches)
end

function _group_indices(cm::ContactMatrix, groups)
    if groups isa Integer
        idx = [Int(groups)]
    elseif groups isa AbstractVector{Bool}
        length(groups) == n_groups(cm) || throw(DimensionMismatch(
            "Boolean target mask length $(length(groups)) does not match groups $(n_groups(cm))"))
        idx = findall(groups)
    elseif groups isa AbstractVector{<:Integer}
        idx = Int.(groups)
    elseif groups isa AbstractString
        idx = [_group_index_from_label(cm, groups)]
    elseif groups isa AbstractVector{<:AbstractString}
        idx = [_group_index_from_label(cm, group) for group in groups]
    else
        throw(ArgumentError("target groups must be indices, labels, or a Boolean mask"))
    end

    isempty(idx) && throw(ArgumentError("target groups must not be empty"))
    allunique(idx) || throw(ArgumentError("target groups must be unique"))
    all(i -> 1 <= i <= n_groups(cm), idx) ||
        throw(ArgumentError("target group indices must be in 1:$(n_groups(cm))"))
    idx
end

function _group_index_from_label(cm::ContactMatrix, label::AbstractString)
    matches = findall(==(String(label)), group_labels(cm))
    length(matches) == 1 || throw(ArgumentError(
        "target label $(repr(label)) must match exactly one group label"))
    only(matches)
end
