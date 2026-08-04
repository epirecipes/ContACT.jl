library(data.table)

normalise_weighted_matrix <- function(target_pop, weighted_matrix) {
  normalised_weighted_matrix <- target_pop * weighted_matrix
  normalised_weighted_matrix <- 0.5 / target_pop *
    (normalised_weighted_matrix + t(normalised_weighted_matrix))
  normalised_weighted_matrix
}

##############################################
# 1. Load ContACT's raw matrix (row = contactee, col = participant)
##############################################
raw <- fread("polymod_uk_age_raw_matrix.csv")
labels <- raw[[1]]
M_contact <- as.matrix(raw[, -1, with = FALSE])
rownames(M_contact) <- labels
colnames(M_contact) <- labels

##############################################
# 2. Reproject in socialmixr's convention (row = participant, col = contactee)
##############################################
M_socialmixr <- t(M_contact)
target_pop <- c(10500.0, 34000.0, 13500.0)   # [0,18), [18,65), 65+
M_reprojected_socialmixr <- normalise_weighted_matrix(target_pop, M_socialmixr)
M_reprojected_contact <- t(M_reprojected_socialmixr)

write.csv(round(M_reprojected_contact, 6), "reprojection_polymod_matrix_R.csv", row.names = TRUE)

cat("\nReciprocity check under target_pop:\n")
for (i in seq_along(labels)) {
  for (j in seq_along(labels)) {
    lhs <- M_reprojected_contact[i, j] * target_pop[j]
    rhs <- M_reprojected_contact[j, i] * target_pop[i]
    cat(sprintf("  (%s,%s): %.6f vs %.6f\n", labels[i], labels[j], lhs, rhs))
  }
}
