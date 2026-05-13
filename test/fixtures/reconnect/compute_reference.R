library(data.table)

part_common <- fread("reconnect_data/reconnect_participant_common.csv")
part_extra <- fread("reconnect_data/reconnect_participant_extra.csv")
sday <- fread("reconnect_data/reconnect_sday.csv")
cnt_common <- fread("reconnect_data/reconnect_contact_common.csv")
cnt_extra <- fread("reconnect_data/reconnect_contact_extra.csv")

part <- merge(part_common, part_extra, by = "part_id")
part <- merge(part, sday, by = "part_id")
cnt <- merge(cnt_common, cnt_extra, by = c("cont_id", "part_id"))

# dayofweek: 0=Monday..6=Sunday (check day_week column from extra)
cat("day_week unique:", paste(sort(unique(part$day_week)), collapse=", "), "\n")
cat("dayofweek unique:", paste(sort(unique(sday$dayofweek)), collapse=", "), "\n")

# Use the day_week column from participant_extra (already weekday/weekend)
day_counts <- part[!is.na(day_week), .N, by = day_week]
cat("\nDay type distribution:\n")
print(day_counts)
day_counts[, sample_prop := N / sum(N)]
day_counts[, target_prop := fifelse(day_week == "weekday", 5/7, 2/7)]
day_counts[, day_weight := target_prop / sample_prop]
cat("\nDay weights:\n")
print(day_counts)
part <- merge(part, day_counts[, .(day_week, day_weight)], by = "day_week")

# How many contacts have known age? 
cat("\nContacts with cnt_age_exact:", sum(!is.na(cnt$cnt_age_exact)), "of", nrow(cnt), "\n")
cat("Contacts with cnt_age_group:", sum(!is.na(cnt$cnt_age_group) & cnt$cnt_age_group != ""), "of", nrow(cnt), "\n")
cat("Contacts with cnt_location:", sum(!is.na(cnt$cnt_location) & cnt$cnt_location != ""), "of", nrow(cnt), "\n")

age_breaks <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, Inf)
age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
                "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")

# Bin ages from exact
part[, p_age_group := cut(part_age_exact, breaks = age_breaks, labels = age_labels, right = FALSE)]
cnt[, c_age_group_exact := cut(cnt_age_exact, breaks = age_breaks, labels = age_labels, right = FALSE)]

##############################################
# 1. Age matrix (individually reported contacts with known exact age, day-weighted)
##############################################
cnt_with_age <- cnt[!is.na(c_age_group_exact)]
contact_counts <- cnt_with_age[, .(n_contacts = .N), by = .(part_id, c_age_group_exact)]
all_combos <- CJ(part_id = unique(part$part_id), c_age_group_exact = factor(age_labels, levels = age_labels))
cc_full <- merge(all_combos, contact_counts, by = c("part_id", "c_age_group_exact"), all.x = TRUE)
cc_full[is.na(n_contacts), n_contacts := 0]
cc_full <- merge(cc_full, part[, .(part_id, p_age_group, day_weight)], by = "part_id")

age_matrix <- cc_full[!is.na(p_age_group), .(
  weighted_mean = sum(n_contacts * day_weight) / sum(day_weight)
), by = .(p_age_group, c_age_group_exact)]

age_wide <- dcast(age_matrix, p_age_group ~ c_age_group_exact, value.var = "weighted_mean")
rn <- age_wide$p_age_group; age_wide[, p_age_group := NULL]
M_age <- as.matrix(age_wide); rownames(M_age) <- rn

cat("\n=== AGE MATRIX (day-weighted, exact age) ===\n")
cat("Row sums:\n"); print(round(rowSums(M_age), 4))
cat("Grand mean total contacts:", round(sum(part$day_weight * sapply(part$part_id, function(pid) sum(cnt_with_age$part_id == pid))) / sum(part$day_weight), 4), "\n")
write.csv(round(M_age, 6), "reconnect_age_matrix_R.csv", row.names = TRUE)

# Participants per age group
part_n <- part[!is.na(p_age_group), .N, by = p_age_group][order(factor(p_age_group, levels = age_labels))]
cat("\nParticipants per age group:\n")
print(part_n)

##############################################
# 2. Ethnicity matrix (complete case, day-weighted)
##############################################
eth_levels <- c("Asian","Black","Mixed","Other","White")
cnt_eth <- cnt[cnt_ethnicity %in% eth_levels]
part_eth <- part[part_ethnicity %in% eth_levels]
eth_counts <- cnt_eth[part_id %in% part_eth$part_id, .(n_contacts = .N), by = .(part_id, cnt_ethnicity)]
eth_combos <- CJ(part_id = unique(part_eth$part_id), cnt_ethnicity = factor(eth_levels, levels = eth_levels))
eth_full <- merge(eth_combos, eth_counts, by = c("part_id", "cnt_ethnicity"), all.x = TRUE)
eth_full[is.na(n_contacts), n_contacts := 0]
eth_full <- merge(eth_full, part_eth[, .(part_id, part_ethnicity, day_weight)], by = "part_id")

eth_matrix <- eth_full[, .(weighted_mean = sum(n_contacts * day_weight) / sum(day_weight)),
                        by = .(part_ethnicity, cnt_ethnicity)]
eth_wide <- dcast(eth_matrix, part_ethnicity ~ cnt_ethnicity, value.var = "weighted_mean")
rn_e <- eth_wide$part_ethnicity; eth_wide[, part_ethnicity := NULL]
M_eth <- as.matrix(eth_wide); rownames(M_eth) <- rn_e
cat("\n=== ETHNICITY MATRIX ===\n"); print(round(M_eth, 4))
write.csv(round(M_eth, 6), "reconnect_ethnicity_matrix_R.csv", row.names = TRUE)

# Ethnicity population proportions (from survey, unweighted)
eth_props <- part_eth[, .N, by = part_ethnicity][order(factor(part_ethnicity, levels = eth_levels))]
eth_props[, prop := N / sum(N)]
cat("\nEthnicity proportions (survey):\n"); print(eth_props)

##############################################
# 3. Setting-specific age matrices
##############################################
for (loc in c("Home","Work","School","Other")) {
  cnt_loc <- cnt[cnt_location == loc & !is.na(c_age_group_exact)]
  loc_counts <- cnt_loc[, .(n_contacts = .N), by = .(part_id, c_age_group_exact)]
  loc_combos <- CJ(part_id = unique(part$part_id), c_age_group_exact = factor(age_labels, levels = age_labels))
  loc_full <- merge(loc_combos, loc_counts, by = c("part_id", "c_age_group_exact"), all.x = TRUE)
  loc_full[is.na(n_contacts), n_contacts := 0]
  loc_full <- merge(loc_full, part[, .(part_id, p_age_group, day_weight)], by = "part_id")
  loc_mat <- loc_full[!is.na(p_age_group), .(wm = sum(n_contacts * day_weight) / sum(day_weight)),
                       by = .(p_age_group, c_age_group_exact)]
  loc_wide <- dcast(loc_mat, p_age_group ~ c_age_group_exact, value.var = "wm")
  rn_l <- loc_wide$p_age_group; loc_wide[, p_age_group := NULL]
  M_loc <- as.matrix(loc_wide); rownames(M_loc) <- rn_l
  write.csv(round(M_loc, 6), paste0("reconnect_age_", tolower(loc), "_matrix_R.csv"), row.names = TRUE)
  cat(paste0("\n", loc, " row sums: ", paste(round(rowSums(M_loc), 3), collapse=", "), "\n"))
}

# Composition check
M_h <- as.matrix(fread("reconnect_age_home_matrix_R.csv", drop = 1))
M_w <- as.matrix(fread("reconnect_age_work_matrix_R.csv", drop = 1))
M_s <- as.matrix(fread("reconnect_age_school_matrix_R.csv", drop = 1))
M_o <- as.matrix(fread("reconnect_age_other_matrix_R.csv", drop = 1))

# Contacts with known location AND known age
cnt_known_loc_age <- cnt[!is.na(cnt_location) & cnt_location != "" & !is.na(c_age_group_exact)]
cnt_all_age <- cnt[!is.na(c_age_group_exact)]
cat("\nContacts with known loc+age:", nrow(cnt_known_loc_age), "vs all with age:", nrow(cnt_all_age), "\n")

##############################################
# 4. Coarse 3x3 
##############################################
coarse_breaks <- c(0, 18, 65, Inf)
coarse_labels <- c("0-17","18-64","65+")
part[, p_age_coarse := cut(part_age_exact, breaks = coarse_breaks, labels = coarse_labels, right = FALSE)]
cnt[, c_age_coarse := cut(cnt_age_exact, breaks = coarse_breaks, labels = coarse_labels, right = FALSE)]

cc3 <- cnt[!is.na(c_age_coarse), .(n = .N), by = .(part_id, c_age_coarse)]
c3_combos <- CJ(part_id = unique(part$part_id), c_age_coarse = factor(coarse_labels, levels = coarse_labels))
c3_full <- merge(c3_combos, cc3, by = c("part_id", "c_age_coarse"), all.x = TRUE)
c3_full[is.na(n), n := 0]
c3_full <- merge(c3_full, part[, .(part_id, p_age_coarse, day_weight)], by = "part_id")

c3_mat <- c3_full[!is.na(p_age_coarse), .(wm = sum(n * day_weight) / sum(day_weight)),
                   by = .(p_age_coarse, c_age_coarse)]
c3_wide <- dcast(c3_mat, p_age_coarse ~ c_age_coarse, value.var = "wm")
rn3 <- c3_wide$p_age_coarse; c3_wide[, p_age_coarse := NULL]
M_coarse <- as.matrix(c3_wide); rownames(M_coarse) <- rn3

cat("\n=== COARSE 3x3 ===\n"); print(round(M_coarse, 4))
write.csv(round(M_coarse, 6), "reconnect_age_coarse_matrix_R.csv", row.names = TRUE)

# Population per coarse group (from survey)
pop_coarse <- part[!is.na(p_age_coarse), .N, by = p_age_coarse][order(factor(p_age_coarse, levels = coarse_labels))]
pop_coarse[, prop := N / sum(N)]
cat("\nCoarse population fractions (survey):\n"); print(pop_coarse)

cat("\nAll done.\n")
