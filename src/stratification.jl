"""
Stratification: replicating contact structure across spatial/demographic strata.

Given a local contact matrix and a coupling (mixing) matrix between strata, the
stratified matrix is block-structured:

    A_strat[(s_to, i), (s_from, j)] = D[s_to, s_from] · A_local[i, j]

which is the Kronecker product A_strat = D ⊗ A_local, bilinear in the two arguments.

This builds a matrix over a product index; it does not stratify an epidemic model.
Model stratification is a separate operation (AlgebraicPetri's `typed_product`); to
pair this matrix with it, align the indices through a `FinFunction`, not by position.
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

Written in `MeanContacts` — replicating an entry across strata leaves mean contacts
per participant unchanged, while the stratum populations are split, so the counts and
per-capita readings do change — and carried to the other representations by the
representation isomorphisms, so that
`reinterpret_units(stratify(cm, D), s) ≈ stratify(reinterpret_units(cm, s), D)` — exact
in real arithmetic, agreeing to about 15 significant digits in floating point.
"""
stratify(cm::ContactMatrix, coupling::AbstractMatrix{<:Real}; kwargs...) =
    _via(MeanContacts(), c -> _stratify_mean(c, coupling; kwargs...), cm)

function _stratify_mean(cm::ContactMatrix, coupling::AbstractMatrix{<:Real};
                  stratum_populations::Union{Nothing, AbstractVector{<:Real}, AbstractMatrix{<:Real}}=nothing,
                  stratum_labels::Union{Nothing, Vector{String}}=nothing)
    M = matrix(cm)
    n_local = n_groups(cm)
    n_strata = size(coupling, 1)
    size(coupling, 1) == size(coupling, 2) || throw(DimensionMismatch(
        "coupling matrix must be square, got $(size(coupling))"))
    _valid_coupling(coupling) ||
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

    strat_partition = stratified_partition(cm.partition, n_strata; stratum_labels)
    strat_pop = _stratum_population(population(cm), n_local, n_strata, stratum_populations)

    ContactMatrix(A_strat, strat_partition, strat_pop, cm.semantics)
end

"""
    stratified_partition(local_partition, n_strata; stratum_labels=nothing)

The partition of a stratified contact matrix: strata × `local_partition`, with strata
labelled `stratum_labels` (default `"S1"`, `"S2"`, …). Index `(s - 1)·n_local + i` is
local group `i` in stratum `s`.

This is the partition `stratify` gives its result, and the domain of the `PartitionMap`
that [`stratify_coarsen`](@ref) takes.
"""
function stratified_partition(local_partition::AbstractPartition, n_strata::Integer;
                              stratum_labels::Union{Nothing, AbstractVector{<:AbstractString}}=nothing)
    labels = stratum_labels === nothing ? ["S$s" for s in 1:n_strata] : String.(stratum_labels)
    length(labels) == n_strata || throw(DimensionMismatch(
        "stratum_labels must have length $n_strata, got $(length(labels))"))
    strata = CategoricalPartition{:stratum,Int}(collect(1:n_strata);
        participant_col=:stratum, contact_col=:stratum, labels=labels)
    strata × local_partition
end

# Replicated population for the single-matrix form. By default, split the local
# population equally across strata; callers can pass explicit expanded populations.
function _stratum_population(pop::AbstractVector, n_local::Integer, n_strata::Integer,
                             stratum_populations)
    n_total = n_local * n_strata
    if stratum_populations === nothing
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
end

"""
    stratify(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix; stratum_labels=nothing)

Heterogeneous stratification: each stratum has its own local contact matrix,
connected by a coupling matrix. The source stratum's matrix is used.

All inputs must share a representation; the result is returned in it, transported
from `MeanContacts` as in the single-matrix method. `stratum_labels` names the strata
(default `"S1"`, `"S2"`, …).
"""
function stratify(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix{<:Real};
                  stratum_labels::Union{Nothing, AbstractVector{<:AbstractString}}=nothing)
    isempty(cms) && throw(ArgumentError("cannot stratify an empty collection"))
    all(cm -> typeof(cm.semantics) == typeof(cms[1].semantics), cms) || throw(ArgumentError(
        "all contact matrices must have the same unit semantics"))
    home = MeanContacts()
    cms[1].semantics isa MeanContacts && return _stratify_mean(cms, coupling; stratum_labels)
    converted = [reinterpret_units(cm, home) for cm in cms]
    reinterpret_units(_stratify_mean(converted, coupling; stratum_labels), cms[1].semantics)
end

function _check_stratify_inputs(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix{<:Real})
    n_strata = length(cms)
    size(coupling) == (n_strata, n_strata) || throw(DimensionMismatch(
        "coupling must be $n_strata × $n_strata, got $(size(coupling))"))
    _valid_coupling(coupling) ||
        throw(ArgumentError("coupling entries must be finite and non-negative"))

    n_local = n_groups(cms[1])
    all(cm -> n_groups(cm) == n_local, cms) || throw(ArgumentError(
        "all contact matrices must have the same number of age groups"))
    all(cm -> same_partition(cm.partition, cms[1].partition), cms) || throw(ArgumentError(
        "all contact matrices must have the same partition"))
    all(cm -> typeof(cm.semantics) == typeof(cms[1].semantics), cms) || throw(ArgumentError(
        "all contact matrices must have the same unit semantics"))
    n_local
end

function _stratify_mean(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix{<:Real};
                        stratum_labels=nothing)
    n_strata = length(cms)
    n_local = _check_stratify_inputs(cms, coupling)

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

    strat_partition = stratified_partition(cms[1].partition, n_strata; stratum_labels)
    pop = vcat([population(cm) for cm in cms]...)
    ContactMatrix(A_strat, strat_partition, pop, cms[1].semantics)
end

# ---------------------------------------------------------------------------
# Blockwise ⊗ then ↓
# ---------------------------------------------------------------------------

"""
    stratify_coarsen(cm::ContactMatrix, coupling, f::PartitionMap; stratum_populations=nothing, stratum_labels=nothing)
    stratify_coarsen(cms::AbstractVector{<:ContactMatrix}, coupling, f::PartitionMap; stratum_labels=nothing)

Stratify across strata connected by a coupling matrix and coarsen along `f` in one step.
The result equals `coarsen(stratify(cm, coupling; …), f)`, or
`coarsen(stratify(cms, coupling), f)`, but the stratified matrix is never formed.

The coarse total contacts are

    T_coarse[I, J] = Σ_{s,r} D[s, r] · Σ_{f(s,i) = I} Σ_{f(r,j) = J} M_r[i, j] · N_r[j],

where `D` is `coupling`, and `M_r`, `N_r` are the local matrix and group populations of
stratum `r`. The work is proportional to the number of nonzero coupling entries times
`n_local²`, and the memory to the size of the result, so many strata can be handled, for
example small areas summed to a national matrix. A `Diagonal` coupling, for strata that
do not mix, is read along its diagonal.

# Arguments
- `cm`, `stratum_populations`, `stratum_labels`: as in the single-matrix `stratify`.
- `cms`: as in the heterogeneous `stratify`, one local matrix per stratum.
- `coupling`: `n_strata × n_strata` non-negative mixing matrix between strata.
- `f`: a `PartitionMap` whose domain is the stratified partition; build it with
  [`stratified_partition`](@ref). For example,
  `PartitionMap(stratified_partition(p, n_strata), p)` sums over strata.

# Returns
A `ContactMatrix` over the codomain of `f`, in the representation of the input.

Throws for the same inputs as `coarsen(stratify(…), f)`, including a zero-population
group that carries contacts along an axis the input's representation rescales. Written
in `MeanContacts` and transported to the other representations; agrees with
`coarsen(stratify(…), f)` to floating-point roundoff.
"""
function stratify_coarsen(cm::ContactMatrix, coupling::AbstractMatrix{<:Real}, f::PartitionMap;
                          stratum_populations::Union{Nothing, AbstractVector{<:Real}, AbstractMatrix{<:Real}}=nothing,
                          stratum_labels::Union{Nothing, AbstractVector{<:AbstractString}}=nothing)
    semantics = cm.semantics
    home_op = function (c)
        n_local = n_groups(c)
        n_strata = size(coupling, 1)
        size(coupling, 1) == size(coupling, 2) || throw(DimensionMismatch(
            "coupling matrix must be square, got $(size(coupling))"))
        _valid_coupling(coupling) ||
            throw(ArgumentError("coupling entries must be finite and non-negative"))
        partition = stratified_partition(c.partition, n_strata; stratum_labels)
        strat_pop = _stratum_population(population(c), n_local, n_strata, stratum_populations)
        pops = [strat_pop[((r - 1) * n_local + 1):(r * n_local)] for r in 1:n_strata]
        blocks = fill(matrix(c), n_strata)
        _stratify_coarsen_mean(blocks, pops, coupling, partition, f, semantics)
    end
    _via(MeanContacts(), home_op, cm)
end

function stratify_coarsen(cms::AbstractVector{<:ContactMatrix}, coupling::AbstractMatrix{<:Real},
                          f::PartitionMap;
                          stratum_labels::Union{Nothing, AbstractVector{<:AbstractString}}=nothing)
    isempty(cms) && throw(ArgumentError("cannot stratify an empty collection"))
    all(cm -> typeof(cm.semantics) == typeof(cms[1].semantics), cms) || throw(ArgumentError(
        "all contact matrices must have the same unit semantics"))
    semantics = cms[1].semantics
    home = MeanContacts()
    converted = semantics isa MeanContacts ? cms : [reinterpret_units(cm, home) for cm in cms]
    _check_stratify_inputs(converted, coupling)
    partition = stratified_partition(converted[1].partition, length(converted); stratum_labels)
    out = _stratify_coarsen_mean([matrix(cm) for cm in converted],
                                 [population(cm) for cm in converted],
                                 coupling, partition, f, semantics)
    semantics isa MeanContacts ? out : reinterpret_units(out, semantics)
end

# Coarse result in MeanContacts from local blocks M_r and stratum populations N_r.
# `semantics` is the caller's representation: the stratified matrix must be expressible
# in it, which `_check_stratified_reinterpretable` checks block by block.
function _stratify_coarsen_mean(blocks, pops, coupling, partition, f::PartitionMap,
                                semantics::UnitSemantics)
    n_strata = length(blocks)
    n_local = size(blocks[1], 1)
    all(x -> isfinite(x) && x >= 0, Iterators.flatten(pops)) ||
        throw(ArgumentError("population entries must be finite and non-negative"))
    _check_stratified_reinterpretable(blocks, pops, coupling, semantics)
    same_partition(partition, f.domain) || throw(ArgumentError(
        "ContactMatrix partition does not match PartitionMap domain"))

    fmap = _assignments(f)
    n_coarse = n_groups(f.codomain)
    index(s, i) = (s - 1) * n_local + i

    pop_coarse = zeros(Float64, n_coarse)
    for r in 1:n_strata, j in 1:n_local
        pop_coarse[fmap[index(r, j)]] += pops[r][j]
    end

    # Total contacts per block: T_r[i, j] = M_r[i, j] · N_r[j].
    totals = [blocks[r] .* permutedims(pops[r]) for r in 1:n_strata]
    T_coarse = zeros(Float64, n_coarse, n_coarse)
    for idx in _nonzero_indices(coupling)
        s, r = Tuple(idx)
        d = coupling[idx]
        Tr = totals[r]
        for j in 1:n_local
            J = fmap[index(r, j)]
            for i in 1:n_local
                T_coarse[fmap[index(s, i)], J] += d * Tr[i, j]
            end
        end
    end

    M_coarse = zeros(Float64, n_coarse, n_coarse)
    for J in 1:n_coarse
        pop_coarse[J] > 0 || continue
        for I in 1:n_coarse
            M_coarse[I, J] = T_coarse[I, J] / pop_coarse[J]
        end
    end
    ContactMatrix(M_coarse, f.codomain, pop_coarse, MeanContacts())
end

_valid_coupling(coupling::AbstractMatrix) = all(x -> isfinite(x) && x >= 0, coupling)
_valid_coupling(coupling::Diagonal) = all(x -> isfinite(x) && x >= 0, coupling.diag)

# Nonzero coupling entries. A Diagonal coupling (strata that do not mix) is read along
# its diagonal.
_nonzero_indices(coupling::AbstractMatrix) = findall(!iszero, coupling)
_nonzero_indices(coupling::Diagonal) =
    [CartesianIndex(k, k) for k in axes(coupling, 1) if !iszero(coupling[k, k])]

# Throws unless the stratified matrix can be converted from MeanContacts to `semantics`.
# As in `reinterpret_units`, a zero-population group may carry no contacts along an axis
# the conversion rescales: rows for a change of the row exponent, columns for a change of
# the column exponent. Row (s, i) of the stratified matrix is D[s, r]·M_r[i, :] over r;
# column (r, j) is D[:, r] ⊗ M_r[:, j].
function _check_stratified_reinterpretable(blocks, pops, coupling, semantics::UnitSemantics)
    a, b = population_exponents(MeanContacts())
    a′, b′ = population_exponents(semantics)
    row_exp, col_exp = a - a′, b - b′
    n_strata = length(blocks)
    n_local = size(blocks[1], 1)
    if row_exp != 0
        for s in 1:n_strata, i in 1:n_local
            pops[s][i] == 0 || continue
            for r in 1:n_strata
                (coupling[s, r] != 0 && !all(iszero, @view blocks[r][i, :])) && throw(ArgumentError(
                    "cannot reinterpret row $((s - 1) * n_local + i): zero population with nonzero contacts"))
            end
        end
    end
    if col_exp != 0
        for r in 1:n_strata, j in 1:n_local
            pops[r][j] == 0 || continue
            (!all(iszero, @view coupling[:, r]) && !all(iszero, @view blocks[r][:, j])) && throw(ArgumentError(
                "cannot reinterpret column $((r - 1) * n_local + j): zero population with nonzero contacts"))
        end
    end
    nothing
end
