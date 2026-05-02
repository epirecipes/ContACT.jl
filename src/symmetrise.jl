"""
Symmetrisation morphism.

Makes a contact matrix symmetric so that c_ij · N_i = c_ji · N_j,
preserving the reciprocity constraint (total contacts from group i to j
equals total contacts from j to i).

Key property (proven in Lean): symmetrise is idempotent.
    symmetrise(symmetrise(M)) = symmetrise(M)
"""

"""
    symmetrise(cm::ContactMatrix)

Symmetrise a contact matrix preserving the reciprocity constraint.

For each pair (i, j), the symmetrised entry is:
    M_sym[i, j] = (M[i, j] · N_j + M[j, i] · N_i) / (2 · N_j)

This ensures c_ij · N_i = c_ji · N_j in the result.

# Returns
A new `ContactMatrix` with a symmetric contact pattern.
"""
function symmetrise(cm::ContactMatrix)
    M = matrix(cm)
    pop = population(cm)
    n = n_groups(cm)

    M_sym = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            # Total contacts from i→j and j→i, averaged
            total_ij = M[i, j] * pop[j] + M[j, i] * pop[i]
            M_sym[i, j] = pop[j] > 0 ? total_ij / (2.0 * pop[j]) : 0.0
        end
    end

    ContactMatrix(M_sym, cm.partition, pop, cm.semantics)
end
