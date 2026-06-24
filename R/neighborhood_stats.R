#' @import data.table
#' @import Giotto
#' @importFrom stats p.adjust sd median
#' @importFrom dplyr group_by summarise mutate pull n arrange desc
#' @importFrom tidyr separate pivot_wider
#' @importFrom igraph graph_from_data_frame simplify triangles V
#' @importFrom utils txtProgressBar setTxtProgressBar
NULL

# Satisfy R CMD check global variables for non-standard evaluation (NSE)
utils::globalVariables(c(
  "unified_cells", "type_int", "N", "round", "unified_int", "orig", "V1",
  "original", "simulations", "enrichm", "p_higher_orig", "p_lower_orig",
  "p.adj_higher", "p.adj_lower", "PI_value", "int_ranking", "from",
  "to_cell_type", "from_cell_type", "row_total", "list_ID", "source",
  "data", "condition", "mean_val", "sd_val", "n_replicates", "mean.diff",
  "pooled_sd", "effect.size", "t1", "t2", "t3", "i.unified_int", "node1",
  "node2", "node3", "type1", "type2", "type3", "pair1", "pair2", "pair3",
  "original_V1", "sim_V1", "p_min", "p_adj", "rankOrder", "labelSubset"
))

#' Calculate Proximity Enrichment using Simulations
#'
#' Computes the enrichment of cell-cell interactions in a spatial network compared to a simulated background.
#'
#' @param gobject A Giotto object.
#' @param spat_unit Spatial unit to use (default: NULL).
#' @param feat_type Feature type to use (default: NULL).
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

  spat_unit <- Giotto::set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
  feat_type <- Giotto::set_default_feat_type(gobject = gobject, spat_unit = spat_unit, feat_type = feat_type)
  adjust_method <- match.arg(adjust_method)

  spatial_network_annot <- Giotto::annotateSpatialNetwork(
    gobject = gobject,
    feat_type = feat_type,
    spat_unit = spat_unit,
    spatial_network_name = spatial_network_name,
    cluster_column = cluster_column
  )

  spatial_network_annot[, unified_cells := paste(sort(c(to, from)), collapse = "--"), by = 1:nrow(spatial_network_annot)]
  spatial_network_annot <- spatial_network_annot[!duplicated(unified_cells)]

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

  table_sim_results <- sample_dt[, .N, by = c("unified_int", "type_int", "round")]

  unique_ints <- unique(table_sim_results[, .(unified_int, type_int)])
  minimum_simulations <- unique_ints[rep(seq_len(nrow(unique_ints)), number_of_simulations), ]
  minimum_simulations[, round := rep(paste0("sim", seq_len(number_of_simulations)), each = nrow(unique_ints))]
  minimum_simulations[, N := 0]

  table_sim_minimum_results <- rbind(table_sim_results, minimum_simulations)
  table_sim_minimum_results[, V1 := sum(N), by = c("unified_int", "type_int", "round")]
  table_sim_results <- unique(table_sim_minimum_results[, .(unified_int, type_int, round, V1)])

  table_sim_results[, orig := "simulations"]
  spatial_network_annot[, round := "original"]

  table_orig_results <- spatial_network_annot[, .N, by = c("unified_int", "type_int", "round")]
  table_orig_results[, orig := "original"]
  data.table::setnames(table_orig_results, old = "N", new = "V1")

  table_results <- rbind(table_orig_results, table_sim_results)

  all_sim_ints <- unique(table_results[orig == "simulations"]$unified_int)
  all_orig_ints := unique(table_results[orig == "original"]$unified_int)

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
#' Takes the raw simulation table from \code{cellProx.sim} and calculates empirical p-values and enrichment scores.
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

  unique_combos <- unique(table_results$unified_int)
  p_high <- numeric(length(unique_combos))
  p_low <- numeric(length(unique_combos))

  for (i in seq_along(unique_combos)) {
    this_combo <- unique_combos[i]
    sub <- table_results[unified_int == this_combo]

    orig_value <- sub[orig == "original"]$V1
    if (length(orig_value) == 0) orig_value <- 0

    sim_values <- sub[orig == "simulations"]$V1

    if (length(sim_values) < number_of_simulations) {
      sim_values <- c(sim_values, rep(0, number_of_simulations - length(sim_values)))
    }

    p_high[i] <- 1 - (sum((orig_value + 1) > (sim_values + 1)) / number_of_simulations)
    p_low[i] <- 1 - (sum((orig_value + 1) < (sim_values + 1)) / number_of_simulations)
  }

  res_pvalue_DT <- data.table::data.table(
    unified_int = unique_combos,
    p_higher_orig = p_high,
    p_lower_orig = p_low
  )

  table_mean_results <- table_results[, .(V1 = mean(V1)), by = c("orig", "unified_int", "type_int")]
  table_mean_results_dc <- data.table::dcast(table_mean_results, type_int + unified_int ~ orig, value.var = "V1")

  if (!"original" %in% names(table_mean_results_dc)) table_mean_results_dc[, original := 0]
  if (!"simulations" %in% names(table_mean_results_dc)) table_mean_results_dc[, simulations := 0]
  table_mean_results_dc[is.na(original), original := 0]

  table_mean_results_dc[, enrichm := log2((original + 1) / (simulations + 1))]
  table_mean_results_dc <- merge(table_mean_results_dc, res_pvalue_DT, by = "unified_int")

  table_mean_results_dc[, p.adj_higher := stats::p.adjust(p_higher_orig, method = adjust_method)]
  table_mean_results_dc[, p.adj_lower := stats::p.adjust(p_lower_orig, method = adjust_method)]

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
#' Aggregates results from multiple samples and calculates subgroup-wide overall enrichment metrics.
#'
#' @param dt A data.table containing combined simulation results from multiple samples (must have `id`, `unified_int`, `round`, `V1` columns).
#' @param number_of_simulations Integer, matches prior permutation depth (default: 200).
#'
#' @return A list from \code{cellProx.calcP}.
#' @export
cellProx.subgroup <- function(dt, number_of_simulations = 200) {

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

  long_dt[, orig := "NA"]
  long_dt[grep("original", round), orig := "original"]
  long_dt[grep("^sim", round), orig := "simulations"]

  long_dt <- long_dt %>%
    tidyr::separate(col = "unified_int", into = c("left", "right"), sep = "--", remove = FALSE) %>%
    data.table::as.data.table()

  long_dt[, type_int := ifelse(left == right, "homo", "hetero")]
  long_dt[, c("left", "right") := NULL]

  out <- cellProx.calcP(long_dt, number_of_simulations = number_of_simulations)
  return(out)
}

#' Run Sample-Level Interaction Matrix Pipeline
#'
#' Iteratively maps microenvironment profiles and groups them by metadata states.
#'
#' @param allInteractions Raw pooled spatial network tables across samples.
#' @param sampleNames Vector of sample identifiers to execute loop processing across.
#'
#' @return A list of per-sample interaction validation summaries.
#' @export
calc.Interactions <- function(allInteractions, sampleNames) {
  perSampleInteraction.enrichm_res <- list()
  for (i in seq_along(sampleNames)) {
    sub_dt <- allInteractions[id %in% sampleNames[i]]
    calc_obj <- cellProx.subgroup(sub_dt)
    perSampleInteraction.enrichm_res[[sampleNames[i]]] <- calc_obj$enrichm_res
  }
  return(perSampleInteraction.enrichm_res)
}

#' Standardize Sample Labels and Meta Conditions
#'
#' Appends study identifiers and experimental group parameters securely to data collections.
#'
#' @param dt List of processed summary data.tables.
#' @param label Character metadata string defining the group tier (e.g., 'stage3', 'recur').
#'
#' @return A flattened unified data.table ready for downstream differential testing.
#' @export
addLabels.fxn <- function(dt, label) {
  for (i in seq_along(dt)) {
    dt[[i]][, list_ID := names(dt)[i]]
  }
  merged_dt <- data.table::rbindlist(dt)
  merged_dt$condition <- label
  return(merged_dt)
}

#' Define Neighborhood Proportions
#'
#' Calculates the proportion of cell types in the neighborhood of every cell.
#'
#' @param annot.list A list of annotated spatial networks.
#'
#' @return A list containing proportion tables, raw counts, and tracking sample IDs.
#' @export
defineNeighborhoods <- function(annot.list) {
  list_ID <- names(annot.list)
  res <- list()
  copymatrix <- list()

  for (i in list_ID) {
    message(paste0("Processing sample: ", i))
    a <- annot.list[[i]][, .(to, from, from_cell_type, to_cell_type)]

    b <- data.table::copy(a)
    data.table::setnames(b, c("to", "from", "from_cell_type", "to_cell_type"),
                         c("from", "to", "to_cell_type", "from_cell_type"))

    combined <- rbind(a, b)
    neighborhood_counts <- combined[, .N, by = .(from, to_cell_type, from_cell_type)]

    neighborhood_matrix <- data.table::dcast(neighborhood_counts,
                                             from + from_cell_type ~ to_cell_type,
                                             value.var = "N",
                                             fill = 0)

    cols_to_remove <- intersect(names(neighborhood_matrix), c("ARTIFACT", "CD4.outlier", "CD8.outlier"))
    if (length(cols_to_remove) > 0) {
      neighborhood_matrix[, (cols_to_remove) := NULL]
    }

    copymatrix[[i]] <- neighborhood_matrix

    working_dt <- data.table::copy(neighborhood_matrix)
    id_cols <- c("from", "from_cell_type")
    count_cols := setdiff(names(working_dt), id_cols)

    working_dt[, row_total := rowSums(.SD), .SDcols = count_cols]
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
    if (chooseNeighbor %in% names(dt)) {
      sub_filter <- dt[from_cell_type == indexCluster & get(chooseNeighbor) > threshold]
      if (nrow(sub_filter) > 0) {
        sub_filter[, list_ID := i]
        filterList[[i]] <- sub_filter
      }
    }
  }
  return(data.table::rbindlist(filterList, fill = TRUE))
}

#' Permutation Test for Neighborhood Features
#'
#' Tests if a specific feature expression profile shifts significantly in response to close spatial neighbors.
#'
#' @param data.tbl A data.table with columns `source` (group label) and `data` (expression value).
#' @param group.0 Label for the first group (e.g., "No.Neighbors").
#' @param group.1 Label for the second group (e.g., "Neighbors").
#' @param num.permutations Integer, number of permutations (default: 10000).
#'
#' @return The calculated two-tailed p-value.
#' @export
permutation.test.neighborhoods <- function(data.tbl, group.0, group.1, num.permutations = 10000) {
  mean_0 <- mean(data.tbl[source %in% group.0, data], na.rm = TRUE)
  mean_1 <- mean(data.tbl[source %in% group.1, data], na.rm = TRUE)
  obs.diff <- mean_0 - mean_1

  message(paste0("Observed difference: ", round(obs.diff, 2)))

  values <- data.tbl$data
  n_total <- length(values)
  n_g0 <- nrow(data.tbl[source %in% group.0])

  permutation.diffs <- numeric(num.permutations)

  for (i in seq_len(num.permutations)) {
    shuffled_idx <- sample(n_total)
    g0_vals <- values[shuffled_idx[1:n_g0]]
    g1_vals := values[shuffled_idx[(n_g0 + 1):n_total]]
    permutation.diffs[i] <- mean(g0_vals, na.rm = TRUE) - mean(g1_vals, na.rm = TRUE)
  }

  p_val <- mean(abs(permutation.diffs) >= abs(obs.diff))
  message(paste0("P value is ", p_val, ". Lowest possible P value is ", 1 / num.permutations))
  return(p_val)
}

#' Calculate Interaction Effect Sizes
#'
#' Compares interaction enrichment between two conditions using a standardized Cohen's d-like score.
#'
#' @param interaction.dt A data.table containing enrichment results (must have `condition`, `enrichm`, `unified_int`).
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
  working.dt[, enrichm := ifelse(enrichm > 0, enrichm, 0)]

  delta_pairwise <- working.dt[, .(
    mean_val = mean(enrichm, na.rm = TRUE),
    sd_val = sd(enrichm, na.rm = TRUE),
    n_replicates = .N
  ), by = c("unified_int", "condition")]

  wide_stats <- data.table::dcast(delta_pairwise, unified_int ~ condition, value.var = c("mean_val", "sd_val"))

  exp_mean <- paste0("mean_val_", condition.exp)
  ctrl_mean <- paste0("mean_val_", condition.ctrl)
  exp_sd <- paste0("sd_val_", condition.exp)
  ctrl_sd := paste0("sd_val_", condition.ctrl)

  wide_stats[, mean.diff := get(exp_mean) - get(ctrl_mean)]
  wide_stats[, pooled_sd := sqrt(0.5 * (get(exp_sd)^2 + get(ctrl_sd)^2))]
  wide_stats[, effect.size := mean.diff / pooled_sd]

  data.table::setorder(wide_stats, -effect.size)

  df <- wide_stats[abs(effect.size) > ef.threshold]
  df <- df[!unified_int %like% "outlier" & !unified_int %like% "ARTIFACT"]

  if (leaveOutTumor) {
    df <- df[!unified_int %like% "Tumor"]
  }

  return(df)
}

#' Run Comprehensive Triplets Evaluation Loop
#'
#' Executes network mapping, pairwise verification, and subsequent higher-order 3-clique testing.
#'
#' @param gobject A Giotto object.
#' @param nsims Number of background random permutations (default: 200).
#' @param spatnet Character label for naming the spatial network.
#' @param clustCol Cell metadata annotation column tracking cluster names.
#' @param threshold Significance cutoff to select anchored pairs.
#' @param spatialk Number of nearest neighbors for KNN construction.
#' @param sample.names String vector of validation coordinates IDs.
#'
#' @return A list of structural significance frameworks evaluated per-sample.
#' @export
triplet.Analysis <- function(gobject, nsims = 200, spatnet = "spatialknn.arbitrary", clustCol = "L2_annotations", threshold = 0.5,
                             spatialk, sample.names) {
  out.calcP <- list()
  for (i in seq_along(sample.names)) {
    message("starting on sample: ", sample.names[i])
    keepCells <- subset(Giotto::pDataDT(gobject), list_ID == sample.names[i])$cell_ID
    gobj.sm <- Giotto::subsetGiotto(gobject = gobject, cell_ids = keepCells)
    gobj.sm <- Giotto::createSpatialKNNnetwork(gobj.sm, k = spatialk, verbose = TRUE, name = spatnet)

    tableResults.knn <- cellProx.sim(gobj.sm, cluster_column = clustCol, adjust_method = "fdr",
                                     number_of_simulations = nsims, spatial_network_name = spatnet)
    pairwise_dt <- cellProx.subgroup(tableResults.knn, number_of_simulations = nsims)[['enrichm_res']]

    sim_res <- cellTriplet.sim(gobj.sm, spatial_network_name = spatnet, cluster_column = clustCol,
                               number_of_simulations = nsims, pairwise_results_dt = pairwise_dt,
                               pairwise_enrich_threshold = threshold)
    out.calcP[[sample.names[i]]] <- cellTriplet.calcP(table_results = sim_res, number_of_simulations = nsims, adjust_method = "fdr")
  }
  return(out.calcP)
}

#' Calculate Spatial Triplet Enrichment Using Permutations
#'
#' Evaluates higher-order structural motifs (3-cliques / triangles) to characterize complex microenvironment niches.
#'
#' @param gobject A Giotto object.
#' @param spat_unit Spatial unit to use (default: NULL).
#' @param feat_type Feature type to use (default: NULL).
#' @param spatial_network_name Name of the spatial network to use.
#' @param cluster_column Column in cell metadata containing cluster annotations.
#' @param number_of_simulations Number of random simulations to perform (default: 200).
#' @param min_cell_abundance Minimum cells required to keep annotation group (default: 0).
#' @param specific_cell_types Character vector to isolate particular lineages.
#' @param pairwise_results_dt Core background summary data frame containing pairwise calculations.
#' @param pairwise_enrich_threshold Metric cutoff defining verified neighborhood structures.
#'
#' @return A data.table monitoring triplet formations across simulation trials.
#' @export
cellTriplet.sim <- function(gobject, spat_unit = NULL, feat_type = NULL, spatial_network_name,
                            cluster_column, number_of_simulations = 200,
                            min_cell_abundance = 0,
                            specific_cell_types = NULL,
                            pairwise_results_dt = NULL,
                            pairwise_enrich_threshold = 1) {

  get_unified_duplet <- function(t1, t2) { paste(pmin(t1, t2), pmax(t1, t2), sep = "--") }
  get_unified_triplet <- function(t1, t2, t3) {
    apply(cbind(t1, t2, t3), 1, function(row) { paste(sort(row), collapse = "--") })
  }

  spat_net <- Giotto::getSpatialNetwork(gobject, name = spatial_network_name)
  spat_net <- spat_net@networkDT
  cell_meta <- Giotto::pDataDT(gobject)

  type_counts <- cell_meta[, .N, by = cluster_column]
  valid_types <- type_counts[N >= min_cell_abundance][[cluster_column]]
  if(!is.null(specific_cell_types)) { valid_types <- intersect(valid_types, specific_cell_types) }

  cell_meta_filtered <- cell_meta[get(cluster_column) %in% valid_types]
  valid_cell_ids <- cell_meta_filtered$cell_ID

  if(!is.null(pairwise_results_dt)) {
    sig_pairs <- pairwise_results_dt[enrichm <= -pairwise_enrich_threshold | enrichm >= pairwise_enrich_threshold]
    valid_enriched_duplets <- unique(sig_pairs$unified_int)
  } else { stop("Pairwise results must be provided.") }

  all.types <- unique(cell_meta_filtered[[cluster_column]])
  dict <- as.data.table(expand.grid(t1 = all.types, t2 = all.types, t3 = all.types, stringsAsFactors = FALSE))
  dict[, unified_int := get_unified_triplet(t1, t2, t3)]
  data.table::setkey(dict, t1, t2, t3)

  g <- igraph::graph_from_data_frame(spat_net[, .(from, to)], directed = FALSE)
  g <- igraph::simplify(g)

  tri_verts <- igraph::triangles(g)
  tri_mat <- matrix(tri_verts, ncol = 3, byrow = TRUE)
  node_names <- igraph::V(g)$name

  triplet_dt <- data.table(
    node1 = node_names[tri_mat[, 1]],
    node2 = node_names[tri_mat[, 2]],
    node3 = node_names[tri_mat[, 3]]
  )

  triplet_dt <- triplet_dt[node1 %in% valid_cell_ids & node2 %in% valid_cell_ids & node3 %in% valid_cell_ids]
  type_map <- setNames(as.character(cell_meta_filtered[[cluster_column]]), cell_meta_filtered$cell_ID)

  triplet_dt[, type1 := type_map[node1]]
  triplet_dt[, type2 := type_map[node2]]
  triplet_dt[, type3 := type_map[node3]]

  triplet_dt[, pair1 := get_unified_duplet(type1, type2)]
  triplet_dt[, pair2 := get_unified_duplet(type2, type3)]
  triplet_dt[, pair3 := get_unified_duplet(type1, type3)]

  triplet_dt <- triplet_dt[pair1 %in% valid_enriched_duplets | pair2 %in% valid_enriched_duplets | pair3 %in% valid_enriched_duplets]
  triplet_dt[dict, unified_int := i.unified_int, on = .(type1 = t1, type2 = t2, type3 = t3)]

  table_orig_results <- triplet_dt[, .N, by = unified_int]
  table_orig_results[, `:=`(round = "original", orig = "original")]
  setnames(table_orig_results, "N", "V1")

  sim_results_list <- vector("list", number_of_simulations)
  node_mat <- as.matrix(triplet_dt[, .(node1, node2, node3)])

  progBar <- utils::txtProgressBar(min = 0, max = number_of_simulations, style = 3)
  for(i in seq_len(number_of_simulations)) {
    utils::setTxtProgressBar(progBar, i)
    shuffled_types <- sample(type_map)
    names(shuffled_types) <- names(type_map)

    sim_temp <- data.table(type1 = shuffled_types[node_mat[, 1]],
                           type2 = shuffled_types[node_mat[, 2]],
                           type3 = shuffled_types[node_mat[, 3]])
    sim_temp[dict, unified_int := i.unified_int, on = .(type1 = t1, type2 = t2, type3 = t3)]
    sim_dt <- sim_temp[, .N, by = unified_int]
    sim_dt[, round := paste0("sim", i)]
    sim_results_list[[i]] <- sim_dt
  }
  close(progBar)

  table_sim_results <- data.table::rbindlist(sim_results_list)
  table_sim_results[, orig := "simulations"]
  setnames(table_sim_results, "N", "V1")

  unique_ints <- unique(c(table_orig_results$unified_int, table_sim_results$unified_int))
  grid_dt <- as.data.table(expand.grid(unified_int = unique_ints, round = paste0("sim", seq_len(number_of_simulations))))

  table_sim_withZeroes <- merge(grid_dt, table_sim_results[, .(unified_int, round, V1)], by = c("unified_int", "round"), all.x = TRUE)
  table_sim_withZeroes[is.na(V1), V1 := 0]
  table_sim_withZeroes[, orig := "simulations"]

  table_results <- rbind(table_orig_results, table_sim_withZeroes, fill = TRUE)
  return(table_results)
}

#' Calculate P-values for Higher-Order Triplet Structures
#'
#' Scores statistical shift profiles for multi-cellular structures from execution traces.
#'
#' @param table_results Data.table tracker generated via \code{cellTriplet.sim}.
#' @param number_of_simulations Total evaluation step counts.
#' @param adjust_method Standard multi-hypothesis testing modifier tool.
#'
#' @return Comprehensive score frameworks cataloging neighborhood composition properties.
#' @export
cellTriplet.calcP <- function(table_results, number_of_simulations = 1000,
                              adjust_method = c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel", "BY")) {
  sel_adjust_method <- match.arg(adjust_method)

  orig_dt <- table_results[orig == "original", .(unified_int, original_V1 = V1)]
  sim_dt  <- table_results[orig == "simulations", .(unified_int, round, sim_V1 = V1)]
  p_comp <- merge(sim_dt, orig_dt, by = "unified_int", all.x = TRUE)
  p_comp[is.na(original_V1), original_V1 := 0]

  res_pvalue_DT <- p_comp[, .(p_min = min((sum(sim_V1 >= original_V1) + 1) / (number_of_simulations + 1),
                                          (sum(sim_V1 <= original_V1) + 1) / (number_of_simulations + 1)))
                          , by = unified_int]

  table_mean_results <- table_results[, .(V1 = mean(V1)), by = .(orig, unified_int)]
  table_mean_results_dc <- data.table::dcast(table_mean_results, unified_int ~ orig, value.var = "V1")

  if (!"original" %in% names(table_mean_results_dc)) table_mean_results_dc[, original := 0]
  if (!"simulations" %in% names(table_mean_results_dc)) table_mean_results_dc[, simulations := 0]

  table_mean_results_dc[, enrichm := log2((original + 1) / (simulations + 1))]
  table_mean_results_dc <- merge(table_mean_results_dc, res_pvalue_DT, by = "unified_int")

  table_mean_results_dc[, p_adj := stats::p.adjust(p_min, method = sel_adjust_method)]
  table_mean_results_dc[, PI_value := -log10(p_adj + (1 / number_of_simulations)) * enrichm]

  data.table::setorder(table_mean_results_dc, -PI_value)
  table_mean_results_dc[, int_ranking := seq_len(.N)]

  return(list(raw_sim_table = table_results, enrichm_res = table_mean_results_dc))
}
