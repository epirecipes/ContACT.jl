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
