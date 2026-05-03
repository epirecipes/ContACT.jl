"""
Symmetrisation morphism.

Makes a contact matrix reciprocal so that c_ij · N_j = c_ji · N_i,
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

This ensures c_ij · N_j = c_ji · N_i in the result.

If a group has zero population, reciprocal finite rates only exist when the
corresponding total contacts are also zero; otherwise an `ArgumentError` is thrown.

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
            if pop[j] == 0
                total_ij == 0 || throw(ArgumentError(
                    "cannot symmetrise pair ($i, $j): group $j has zero population but nonzero total contacts"))
                M_sym[i, j] = 0.0
            else
                M_sym[i, j] = total_ij / (2.0 * pop[j])
            end
        end
    end

    ContactMatrix(M_sym, cm.partition, pop, cm.semantics)
end
