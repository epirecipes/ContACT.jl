"""
Additive composition of contact matrices.

In the categorical framework, setting-specific contact matrices (home, work,
school, other) combine additively: the total contact matrix is the sum of
setting-specific matrices. This forms a commutative monoid on ContactMatrix
objects sharing the same partition and semantics.

Key properties (proven in Lean):
- Associativity: (A ⊕ B) ⊕ C = A ⊕ (B ⊕ C)
- Commutativity: A ⊕ B = B ⊕ A
- Identity: A ⊕ 0 = A (zero matrix is the identity)
"""

"""
    compose_matrices(a::ContactMatrix, b::ContactMatrix)

Additively compose two contact matrices. Both must share the same
age partition, population, and unit semantics.

Returns a new `ContactMatrix` whose matrix is the elementwise sum.
"""
function compose_matrices(a::ContactMatrix{T1, S}, b::ContactMatrix{T2, S}) where {T1, T2, S}
    a.partition.limits == b.partition.limits || throw(ArgumentError(
        "cannot compose matrices with different age partitions"))
    isapprox(a.population, b.population; rtol=1e-10, atol=0) || throw(ArgumentError(
        "cannot compose matrices with different population vectors"))

    T = promote_type(T1, T2)
    M = Matrix{T}(matrix(a)) + Matrix{T}(matrix(b))
    ContactMatrix(M, a.partition, Vector{T}(a.population), a.semantics)
end

"""
    compose_matrices(matrices::AbstractVector{<:ContactMatrix})

Compose multiple contact matrices additively.
"""
function compose_matrices(matrices::AbstractVector{<:ContactMatrix})
    isempty(matrices) && throw(ArgumentError("cannot compose empty collection"))
    result = matrices[1]
    for i in 2:length(matrices)
        result = compose_matrices(result, matrices[i])
    end
    return result
end

"""
    compose_matrices(; kwargs...)

Compose named contact matrices. Useful for documenting the contribution
of each setting.

# Example
```julia
total = compose_matrices(home=cm_home, work=cm_work, school=cm_school, other=cm_other)
```
"""
function compose_matrices(; kwargs...)
    matrices = collect(values(kwargs))
    compose_matrices(matrices)
end
