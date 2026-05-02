"""
Stratification: replicating contact structure across spatial/demographic strata.

Categorically, stratification is a pullback in a slice category. Given a local
contact matrix and a coupling (mixing) matrix between strata, the stratified
matrix is a block-structured matrix.

    A_strat[(s_to, i), (s_from, j)] = D[s_to, s_from] · A_local[i, j]

This is equivalent to the Kronecker product: A_strat = D ⊗ A_local.
"""

"""
    stratify(cm::ContactMatrix, coupling::AbstractMatrix)

Stratify a contact matrix across strata connected by a coupling matrix.

The resulting block-structured matrix has dimension (n_strata × n_ages).
Every stratum shares the same local contact pattern.

# Arguments
- `cm`: local contact matrix (shared across all strata)
- `coupling`: n_strata × n_strata mixing matrix between strata

# Returns
A `ContactMatrix` with expanded age partition (strata × age groups).
"""
function stratify(cm::ContactMatrix, coupling::AbstractMatrix{<:Real})
    M = matrix(cm)
    n_ages = n_groups(cm)
    n_strata = size(coupling, 1)
    size(coupling, 1) == size(coupling, 2) || throw(DimensionMismatch(
        "coupling matrix must be square, got $(size(coupling))"))

    # Kronecker product: coupling ⊗ local
    n_total = n_strata * n_ages
    T = promote_type(eltype(M), eltype(coupling))
    A_strat = zeros(T, n_total, n_total)

    for s_to in 1:n_strata
        for s_from in 1:n_strata
            rows = ((s_to - 1) * n_ages + 1):(s_to * n_ages)
            cols = ((s_from - 1) * n_ages + 1):(s_from * n_ages)
            A_strat[rows, cols] .= coupling[s_to, s_from] .* M
        end
    end

    # Expanded partition labels
    local_labels = age_labels(cm)
    strat_labels = String[]
    for s in 1:n_strata
        for l in local_labels
            push!(strat_labels, "S$(s):$(l)")
        end
    end
    strat_limits = collect(1.0:n_total)  # synthetic limits for stratified partition
    strat_partition = AgePartition(strat_limits; labels=strat_labels)

    # Replicated population
    pop = population(cm)
    strat_pop = repeat(pop, n_strata) ./ n_strata

    ContactMatrix(A_strat, strat_partition, strat_pop, cm.semantics)
end

"""
    stratify(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix)

Heterogeneous stratification: each stratum has its own local contact matrix,
connected by a coupling matrix. The source stratum's matrix is used.
"""
function stratify(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix{<:Real})
    n_strata = length(cms)
    size(coupling) == (n_strata, n_strata) || throw(DimensionMismatch(
        "coupling must be $n_strata × $n_strata, got $(size(coupling))"))

    n_ages = n_groups(cms[1])
    all(cm -> n_groups(cm) == n_ages, cms) || throw(ArgumentError(
        "all contact matrices must have the same number of age groups"))

    n_total = n_strata * n_ages
    T = promote_type(eltype(matrix(cms[1])), eltype(coupling))
    A_strat = zeros(T, n_total, n_total)

    for s_to in 1:n_strata
        for s_from in 1:n_strata
            rows = ((s_to - 1) * n_ages + 1):(s_to * n_ages)
            cols = ((s_from - 1) * n_ages + 1):(s_from * n_ages)
            A_strat[rows, cols] .= coupling[s_to, s_from] .* matrix(cms[s_from])
        end
    end

    local_labels = age_labels(cms[1])
    strat_labels = String[]
    for s in 1:n_strata
        for l in local_labels
            push!(strat_labels, "S$(s):$(l)")
        end
    end
    strat_limits = collect(1.0:n_total)
    strat_partition = AgePartition(strat_limits; labels=strat_labels)

    pop = vcat([population(cm) for cm in cms]...)
    ContactMatrix(A_strat, strat_partition, pop, cms[1].semantics)
end
