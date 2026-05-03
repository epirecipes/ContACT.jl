"""
Stratification: replicating contact structure across spatial/demographic strata.

Categorically, stratification is a pullback in a slice category. Given a local
contact matrix and a coupling (mixing) matrix between strata, the stratified
matrix is a block-structured matrix.

    A_strat[(s_to, i), (s_from, j)] = D[s_to, s_from] · A_local[i, j]

This is equivalent to the Kronecker product: A_strat = D ⊗ A_local.
"""

"""
    stratify(cm::ContactMatrix, coupling::AbstractMatrix; stratum_populations=nothing, stratum_labels=nothing)

Stratify a contact matrix across strata connected by a coupling matrix.

The resulting block-structured matrix has dimension (n_strata × n_ages).
Every stratum shares the same local contact pattern.

# Arguments
- `cm`: local contact matrix (shared across all strata)
- `coupling`: n_strata × n_strata mixing matrix between strata
- `stratum_populations`: optional expanded population vector, or an
  n_groups × n_strata matrix whose columns are stratum-specific group populations.
  If omitted, the local population is split equally across strata.
- `stratum_labels`: optional labels for the new stratum partition.

# Returns
A `ContactMatrix` with expanded age partition (strata × age groups).
"""
function stratify(cm::ContactMatrix, coupling::AbstractMatrix{<:Real};
                  stratum_populations::Union{Nothing, AbstractVector{<:Real}, AbstractMatrix{<:Real}}=nothing,
                  stratum_labels::Union{Nothing, Vector{String}}=nothing)
    M = matrix(cm)
    n_local = n_groups(cm)
    n_strata = size(coupling, 1)
    size(coupling, 1) == size(coupling, 2) || throw(DimensionMismatch(
        "coupling matrix must be square, got $(size(coupling))"))
    all(x -> isfinite(x) && x >= 0, coupling) ||
        throw(ArgumentError("coupling entries must be finite and non-negative"))

    # Kronecker product: coupling ⊗ local
    n_total = n_strata * n_local
    T = promote_type(eltype(M), eltype(coupling))
    A_strat = zeros(T, n_total, n_total)

    for s_to in 1:n_strata
        for s_from in 1:n_strata
            rows = ((s_to - 1) * n_local + 1):(s_to * n_local)
            cols = ((s_from - 1) * n_local + 1):(s_from * n_local)
            A_strat[rows, cols] .= coupling[s_to, s_from] .* M
        end
    end

    labels = stratum_labels === nothing ? ["S$s" for s in 1:n_strata] : stratum_labels
    length(labels) == n_strata || throw(DimensionMismatch(
        "stratum_labels must have length $n_strata, got $(length(labels))"))
    strata = CategoricalPartition{:stratum,Int}(collect(1:n_strata);
        participant_col=:stratum, contact_col=:stratum, labels=labels)
    strat_partition = strata × cm.partition

    # Replicated population. By default, split the local population equally
    # across strata; callers can pass explicit expanded populations when known.
    pop = population(cm)
    strat_pop = if stratum_populations === nothing
        repeat(pop, n_strata) ./ n_strata
    elseif stratum_populations isa AbstractMatrix
        size(stratum_populations) == (n_local, n_strata) || throw(DimensionMismatch(
            "stratum_populations matrix must be $n_local × $n_strata, got $(size(stratum_populations))"))
        vec(Float64.(stratum_populations))
    else
        length(stratum_populations) == n_total || throw(DimensionMismatch(
            "stratum_populations vector must have length $n_total, got $(length(stratum_populations))"))
        Float64.(stratum_populations)
    end

    ContactMatrix(A_strat, strat_partition, strat_pop, cm.semantics)
end

"""
    stratify(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix)

Heterogeneous stratification: each stratum has its own local contact matrix,
connected by a coupling matrix. The source stratum's matrix is used.
"""
function stratify(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix{<:Real})
    n_strata = length(cms)
    isempty(cms) && throw(ArgumentError("cannot stratify an empty collection"))
    size(coupling) == (n_strata, n_strata) || throw(DimensionMismatch(
        "coupling must be $n_strata × $n_strata, got $(size(coupling))"))
    all(x -> isfinite(x) && x >= 0, coupling) ||
        throw(ArgumentError("coupling entries must be finite and non-negative"))

    n_local = n_groups(cms[1])
    all(cm -> n_groups(cm) == n_local, cms) || throw(ArgumentError(
        "all contact matrices must have the same number of age groups"))
    all(cm -> same_partition(cm.partition, cms[1].partition), cms) || throw(ArgumentError(
        "all contact matrices must have the same partition"))
    all(cm -> typeof(cm.semantics) == typeof(cms[1].semantics), cms) || throw(ArgumentError(
        "all contact matrices must have the same unit semantics"))

    n_total = n_strata * n_local
    T = promote_type((eltype(matrix(cm)) for cm in cms)..., eltype(coupling))
    A_strat = zeros(T, n_total, n_total)

    for s_to in 1:n_strata
        for s_from in 1:n_strata
            rows = ((s_to - 1) * n_local + 1):(s_to * n_local)
            cols = ((s_from - 1) * n_local + 1):(s_from * n_local)
            A_strat[rows, cols] .= coupling[s_to, s_from] .* matrix(cms[s_from])
        end
    end

    strata = CategoricalPartition{:stratum,Int}(collect(1:n_strata);
        participant_col=:stratum, contact_col=:stratum,
        labels=["S$s" for s in 1:n_strata])
    strat_partition = strata × cms[1].partition

    pop = vcat([population(cm) for cm in cms]...)
    ContactMatrix(A_strat, strat_partition, pop, cms[1].semantics)
end
