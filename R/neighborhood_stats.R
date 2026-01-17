#' @import data.table
#' @import Giotto
#' @importFrom stats p.adjust sd median
#' @importFrom dplyr group_by summarise mutate pull n arrange desc
#' @importFrom tidyr separate pivot_wider
NULL

#' Calculate Proximity Enrichment using Simulations
#'
#' Computes the enrichment of cell-cell interactions in a spatial network compared to a simulated background.
#'
#' @param gobject A Giotto object.
#' @param spat_unit Spatial unit to use (default: NULL, uses default).
#' @param feat_type Feature type to use (default: NULL, uses default).
#' @param spatial_network_name Name of the spatial network to use (e.g., "spatialknn.k5").
#' @param cluster_column Column in cell metadata containing cluster/cell-type annotations.
#' @param number_of_simulations Number of random simulations to perform (default: 1000).
#' @param adjust_method Method for p-value adjustment (default: "fdr").
#' @param set_seed Logical, whether to set a seed for reproducibility (default: TRUE).
#' @param seed_number Seed number (default: 1234).
#'
#' @return A data.table containing the counts of interactions for original and simulated networks.
#' @export
cellProx.sim <- function(gobject,
                         spat_unit = NULL,
                         feat_type = NULL,
                         spatial_network_name = "spatialknn.k5",
                         cluster_column,
                         number_of_simulations = 1000,
                         adjust_method = c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel", "BY"),
                         set_seed = TRUE,
                         seed_number = 1234) {
  
  # Set defaults
  spat_unit <- Giotto::set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
  feat_type <- Giotto::set_default_feat_type(gobject = gobject, spat_unit = spat_unit, feat_type = feat_type)
  adjust_method <- match.arg(adjust_method)
  
  # Annotate the spatial network
  # Note: This requires the Giotto object to already have the network created
  spatial_network_annot <- Giotto::annotateSpatialNetwork(
    gobject = gobject,
    feat_type = feat_type,
    spat_unit = spat_unit,
    spatial_network_name = spatial_network_name,
    cluster_column = cluster_column
  )
  
  # Helper for unified cells column (using data.table logic)
  # Note: mimicking GiottoUtils:::dt_sort_combine_two_columns behavior if not exported
  spatial_network_annot[, unified_cells := paste(sort(c(to, from)), collapse = "--"), by = 1:nrow(spatial_network_annot)]
  spatial_network_annot <- spatial_network_annot[!duplicated(unified_cells)]
  
  # Create simulations (using Giotto internal if necessary, or exposed function)
  # WARNING: Giotto:::make_simulated_network is internal. Ensure Giotto version compatibility.
  sample_dt <- Giotto:::make_simulated_network(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type,
    spatial_network_name = spatial_network_name,
    cluster_column = cluster_column,
    number_of_simulations = number_of_simulations,
    set_seed = set_seed,
    seed_number = seed_number
  )
  
  # Combine original and simulated network counts
  table_sim_results <- sample_dt[, .N, by = c("unified_int", "type_int", "round")]
  
  # Fill in zeros for missing simulations
  unique_ints <- unique(table_sim_results[, .(unified_int, type_int)])
  minimum_simulations <- unique_ints[rep(seq_len(nrow(unique_ints)), number_of_simulations), ]
  minimum_simulations[, round := rep(paste0("sim", seq_len(number_of_simulations)), each = nrow(unique_ints))]
  minimum_simulations[, N := 0]
  
  table_sim_minimum_results <- rbind(table_sim_results, minimum_simulations)
  table_sim_minimum_results[, V1 := sum(N), by = c("unified_int", "type_int", "round")]
  table_sim_results <- unique(table_sim_minimum_results[, .(unified_int, type_int, round, V1)])
  
  # Combine with original data
  table_sim_results[, orig := "simulations"]
  spatial_network_annot[, round := "original"]
  
  table_orig_results <- spatial_network_annot[, .N, by = c("unified_int", "type_int", "round")]
  table_orig_results[, orig := "original"]
  data.table::setnames(table_orig_results, old = "N", new = "V1")
  
  table_results <- rbind(table_orig_results, table_sim_results)
  
  # Handle missing combinations in original vs simulations
  all_sim_ints <- unique(table_results[orig == "simulations"]$unified_int)
  all_orig_ints <- unique(table_results[orig == "original"]$unified_int)
  
  missing_in_orig <- setdiff(all_sim_ints, all_orig_ints)
  missing_in_sim <- setdiff(all_orig_ints, all_sim_ints)
  
  if (length(missing_in_orig) > 0) {
    create_missing_orig <- table_results[unified_int %in% missing_in_orig]
    create_missing_orig <- unique(create_missing_orig[, c("orig", "V1") := list("original", 0)])
    create_missing_orig[, round := "original"]
    table_results <- rbind(table_results, create_missing_orig)
  }
  
  if (length(missing_in_sim) > 0) {
    create_missing_sim <- table_results[unified_int %in% missing_in_sim]
    create_missing_sim <- unique(create_missing_sim[, c("orig", "V1") := list("simulations", 0)])
    table_results <- rbind(table_results, create_missing_sim)
  }
  
  return(table_results)
}


#' Calculate P-values for Proximity Enrichment
#'
#' Takes the raw simulation table from \code{cellProx.sim} and calculates p-values and enrichment scores.
#'
#' @param table_results Data.table output from \code{cellProx.sim}.
#' @param number_of_simulations Integer, number of simulations used.
#' @param adjust_method P-value adjustment method.
#' @param set_seed Logical.
#' @param seed_number Integer.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{raw_sim_table}: The input table.
#'   \item \code{enrichm_res}: Data.table with enrichment scores (log2FC) and adjusted p-values.
#' }
#' @export
cellProx.calcP <- function(table_results,
                           number_of_simulations = 1000,
                           adjust_method = "fdr",
                           set_seed = TRUE,
                           seed_number = 1234) {
  
  adjust_method <- match.arg(adjust_method, choices = c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel", "BY"))
  
  # Calculate p-values per interaction combo
  unique_combos <- unique(table_results$unified_int)
  p_high <- numeric(length(unique_combos))
  p_low <- numeric(length(unique_combos))
  
  # Pre-split data for speed
  # Using .SD for subsets might be slower in loop, iterating indices
  for (i in seq_along(unique_combos)) {
    this_combo <- unique_combos[i]
    sub <- table_results[unified_int == this_combo]
    
    orig_value <- sub[orig == "original"]$V1
    if (length(orig_value) == 0) orig_value <- 0
    
    sim_values <- sub[orig == "simulations"]$V1
    
    # Ensure simulation vector length is correct (pad with 0s if missing)
    if (length(sim_values) < number_of_simulations) {
      sim_values <- c(sim_values, rep(0, number_of_simulations - length(sim_values)))
    }
    
    # P-value calculation (empirical)
    p_high[i] <- 1 - (sum((orig_value + 1) > (sim_values + 1)) / number_of_simulations)
    p_low[i] <- 1 - (sum((orig_value + 1) < (sim_values + 1)) / number_of_simulations)
  }
  
  res_pvalue_DT <- data.table::data.table(
    unified_int = unique_combos,
    p_higher_orig = p_high,
    p_lower_orig = p_low
  )
  
  # Calculate mean enrichment
  table_mean_results <- table_results[, .(V1 = mean(V1)), by = c("orig", "unified_int", "type_int")]
  table_mean_results_dc <- data.table::dcast(table_mean_results, type_int + unified_int ~ orig, value.var = "V1")
  
  # Handle NAs
  if (!"original" %in% names(table_mean_results_dc)) table_mean_results_dc[, original := 0]
  if (!"simulations" %in% names(table_mean_results_dc)) table_mean_results_dc[, simulations := 0]
  table_mean_results_dc[is.na(original), original := 0]
  
  # Log2 Enrichment
  table_mean_results_dc[, enrichm := log2((original + 1) / (simulations + 1))]
  
  # Merge with P-values
  table_mean_results_dc <- merge(table_mean_results_dc, res_pvalue_DT, by = "unified_int")
  
  # Adjust P-values
  table_mean_results_dc[, p.adj_higher := stats::p.adjust(p_higher_orig, method = adjust_method)]
  table_mean_results_dc[, p.adj_lower := stats::p.adjust(p_lower_orig, method = adjust_method)]
  
  # Calculate PI value (ranking score)
  table_mean_results_dc[, PI_value := ifelse(p.adj_higher <= p.adj_lower,
                                             -log10(p.adj_higher + (1 / number_of_simulations)) * enrichm,
                                             -log10(p.adj_lower + (1 / number_of_simulations)) * enrichm
  )]
  
  data.table::setorder(table_mean_results_dc, -PI_value)
  table_mean_results_dc[, int_ranking := seq_len(.N)]
  
  return(list(
    raw_sim_table = table_results,
    enrichm_res = table_mean_results_dc
  ))
}


#' Analyze Proximity Subgroups
#'
#' Aggregates results from multiple samples (rbind-ed tableResults) and calculates overall enrichment.
#'
#' @param dt A data.table containing combined simulation results from multiple samples (must have `id`, `unified_int`, `round`, `V1` columns).
#' @param number_of_simulations Integer.
#'
#' @return A list from \code{cellProx.calcP}.
#' @export
cellProx.subgroup <- function(dt, number_of_simulations = 200) {
  
  # Aggregate values across the subgroup
  # We use dcast to sum/mean across the group, then melt back
  wide_dt <- data.table::dcast(dt,
                               unified_int ~ round,
                               value.var = "V1",
                               fun.aggregate = mean,
                               na.rm = TRUE,
                               fill = 0)
  
  long_dt <- data.table::melt(wide_dt,
                              id.vars = "unified_int",
                              variable.name = "round",
                              value.name = "V1",
                              variable.factor = FALSE)
  
  # Recover simulation vs original tags
  long_dt[, orig := "NA"]
  long_dt[grep("original", round), orig := "original"]
  long_dt[grep("^sim", round), orig := "simulations"]
  
  # Parse interaction type (Homo/Hetero)
  long_dt <- long_dt %>%
    tidyr::separate(col = "unified_int", into = c("left", "right"), sep = "--", remove = FALSE) %>%
    data.table::as.data.table()
  
  long_dt[, type_int := ifelse(left == right, "homo", "hetero")]
  long_dt[, c("left", "right") := NULL]
  
  # Calculate P-values for the aggregated group
  out <- cellProx.calcP(long_dt, number_of_simulations = number_of_simulations)
  return(out)
}


#' Define Neighborhood Proportions
#'
#' Calculates the proportion of cell types in the neighborhood of every cell.
#'
#' @param annot.list A list of annotated spatial networks (output of `annotateSpatialNetwork` loop).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{prop.table}: List of data.tables (one per sample) with neighbor proportions.
#'   \item \code{raw.table}: List of raw count matrices.
#'   \item \code{list_ID}: Vector of sample IDs.
#' }
#' @export
defineNeighborhoods <- function(annot.list) {
  list_ID <- names(annot.list)
  res <- list()
  copymatrix <- list()
  
  for (i in list_ID) {
    message(paste0("Processing sample: ", i))
    
    # Network is undirected, so we must count A->B and B->A to get full neighborhood
    a <- annot.list[[i]][, c("to", "from", "from_cell_type", "to_cell_type"), with = FALSE]
    
    # Swap to create the reverse edges
    b <- data.table::copy(a)
    data.table::setnames(b, c("to", "from", "from_cell_type", "to_cell_type"),
                         c("from", "to", "to_cell_type", "from_cell_type"))
    
    combined <- rbind(a, b)
    
    # Count neighbors per cell
    neighborhood_counts <- combined[, .N, by = .(from, to_cell_type, from_cell_type)]
    
    # Cast to wide format (rows = cells, cols = neighbor types)
    neighborhood_matrix <- data.table::dcast(neighborhood_counts,
                                             from + from_cell_type ~ to_cell_type,
                                             value.var = "N",
                                             fill = 0)
    
    # Remove artifacts/outliers if present
    cols_to_remove <- intersect(names(neighborhood_matrix), c("ARTIFACT", "CD4.outlier", "CD8.outlier"))
    if (length(cols_to_remove) > 0) {
      neighborhood_matrix[, (cols_to_remove) := NULL]
    }
    
    copymatrix[[i]] <- neighborhood_matrix
    
    # Calculate proportions
    working_dt <- data.table::copy(neighborhood_matrix)
    id_cols <- c("from", "from_cell_type")
    count_cols <- setdiff(names(working_dt), id_cols)
    
    # Row sums
    working_dt[, row_total := rowSums(.SD), .SDcols = count_cols]
    
    # Divide by total
    working_dt[, (count_cols) := lapply(.SD, function(x) x / row_total), .SDcols = count_cols]
    
    res[[i]] <- working_dt
  }
  
  return(list(prop.table = res, raw.table = copymatrix, list_ID = list_ID))
}


#' Filter Neighborhoods
#'
#' Filters cells based on the proportion of a specific neighbor type.
#'
#' @param res.list Output from \code{defineNeighborhoods}.
#' @param indexCluster The cell type of the "center" cell to filter for.
#' @param chooseNeighbor The cell type of the neighbor to check proportions against.
#' @param threshold The minimum proportion required.
#'
#' @return A combined data.table of cells meeting the criteria.
#' @export
exampleNeighborhoods.filter <- function(res.list, indexCluster, chooseNeighbor, threshold) {
  filterList <- list()
  ids <- res.list[["list_ID"]]
  prop_tables <- res.list[["prop.table"]]
  
  for (i in ids) {
    dt <- prop_tables[[i]]
    # Check if neighbor column exists
    if (chooseNeighbor %in% names(dt)) {
      filterList[[i]] <- dt[from_cell_type == indexCluster & get(chooseNeighbor) > threshold]
      if (nrow(filterList[[i]]) > 0) filterList[[i]][, list_ID := i]
    }
  }
  return(data.table::rbindlist(filterList, fill = TRUE))
}


#' Permutation Test for Neighborhood Features
#'
#' Tests if a specific feature (gene/protein) is differentially expressed in cells that have neighbors of a certain type vs those that don't.
#'
#' @param data.tbl A data.table with columns `source` (group label) and `data` (expression value).
#' @param group.0 Label for the first group (e.g., "No.Neighbors").
#' @param group.1 Label for the second group (e.g., "Neighbors").
#' @param num.permutations Integer, number of permutations.
#'
#' @return The calculated p-value.
#' @export
permutation.test.neighborhoods <- function(data.tbl, group.0, group.1, num.permutations = 10000) {
  
  # Observed difference
  mean_0 <- mean(data.tbl[source %in% group.0, data], na.rm = TRUE)
  mean_1 <- mean(data.tbl[source %in% group.1, data], na.rm = TRUE)
  obs.diff <- mean_0 - mean_1
  
  message(paste0("Observed difference: ", round(obs.diff, 2)))
  
  # Vectorized permutation test
  # Instead of a slow loop, we can create a matrix of shuffles if memory allows,
  # or use a simplified loop. Given 10k perms, a loop is acceptable in C++ but slow in R.
  # We will stick to a clean loop for readability but optimize calculation.
  
  values <- data.tbl$data
  n_total <- length(values)
  n_g0 <- nrow(data.tbl[source %in% group.0])
  
  permutation.diffs <- numeric(num.permutations)
  
  for (i in seq_len(num.permutations)) {
    # Shuffle indices
    shuffled_idx <- sample(n_total)
    # First n_g0 are group 0, rest are group 1
    g0_vals <- values[shuffled_idx[1:n_g0]]
    g1_vals <- values[shuffled_idx[(n_g0 + 1):n_total]]
    permutation.diffs[i] <- mean(g0_vals, na.rm = TRUE) - mean(g1_vals, na.rm = TRUE)
  }
  
  # Two-tailed p-value
  p_val <- mean(abs(permutation.diffs) >= abs(obs.diff))
  
  message(paste0("P value is ", p_val, ". Lowest possible P value is ", 1 / num.permutations))
  return(p_val)
}


#' Calculate Interaction Effect Sizes
#'
#' Compares interaction enrichment between two conditions (e.g., Stage 3 vs Stage 1).
#'
#' @param interaction.dt A data.table containing enrichment results (output of `cellProx.subgroup` for multiple samples). Must have `condition`, `enrichm`, `unified_int`.
#' @param condition.exp Experimental condition string (e.g. "stage3").
#' @param condition.ctrl Control condition string (e.g. "stage1").
#' @param ef.threshold Effect size threshold for filtering.
#' @param leaveOutTumor Logical, whether to exclude tumor interactions.
#'
#' @return A data.table of interactions sorted by effect size.
#' @export
calculateInteractionEffectSizes <- function(interaction.dt,
                                            condition.exp,
                                            condition.ctrl,
                                            ef.threshold = 1,
                                            leaveOutTumor = TRUE) {
  
  working.dt <- data.table::copy(interaction.dt)
  
  # Clip negative enrichment to 0 (per original logic)
  working.dt[, enrichm := ifelse(enrichm > 0, enrichm, 0)]
  
  # Calculate summary stats per condition
  delta_pairwise <- working.dt[, .(
    mean_val = mean(enrichm, na.rm = TRUE),
    sd_val = sd(enrichm, na.rm = TRUE),
    n_replicates = .N
  ), by = c("unified_int", "condition")]
  
  # Pivot to wide to compare conditions
  wide_stats <- data.table::dcast(delta_pairwise, unified_int ~ condition, value.var = c("mean_val", "sd_val"))
  
  exp_mean <- paste0("mean_val_", condition.exp)
  ctrl_mean <- paste0("mean_val_", condition.ctrl)
  exp_sd <- paste0("sd_val_", condition.exp)
  ctrl_sd <- paste0("sd_val_", condition.ctrl)
  
  # Calculate Effect Size (Cohen's d-like)
  # (mean_exp - mean_ctrl) / sqrt(0.5 * (sd_exp^2 + sd_ctrl^2))
  wide_stats[, mean.diff := get(exp_mean) - get(ctrl_mean)]
  wide_stats[, pooled_sd := sqrt(0.5 * (get(exp_sd)^2 + get(ctrl_sd)^2))]
  wide_stats[, effect.size := mean.diff / pooled_sd]
  
  data.table::setorder(wide_stats, -effect.size)
  
  # Filtering
  df <- wide_stats[abs(effect.size) > ef.threshold]
  df <- df[!unified_int %like% "outlier" & !unified_int %like% "ARTIFACT"]
  
  if (leaveOutTumor) {
    df <- df[!unified_int %like% "Tumor"]
  }
  
  return(df)
}


