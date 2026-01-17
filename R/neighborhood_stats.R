
cellProx.sim <- function(gobject, spat_unit = NULL, feat_type = NULL, spatial_network_name="spatialknn.k5", 
                         #spatial_network_name = "Delaunay_network",
                         cluster_column, number_of_simulations = 1000, 
                         adjust_method = c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel", "BY"), 
                         set_seed = TRUE, seed_number = 1234) 
{
  # Set feat_type and spat_unit
  spat_unit <- set_default_spat_unit(
    gobject = gobject,
    spat_unit = spat_unit
  )
  feat_type <- set_default_feat_type(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type
  )
  
  # p.adj test
  sel_adjust_method <- match.arg(adjust_method, choices = c(
    "none", "fdr", "bonferroni", "BH",
    "holm", "hochberg", "hommel",
    "BY"
  ))
  
  spatial_network_annot <- annotateSpatialNetwork(
    gobject = gobject,
    feat_type = feat_type,
    spat_unit = spat_unit,
    spatial_network_name = spatial_network_name,
    cluster_column = cluster_column
  )
  
  # remove double edges between same cells #
  # a simplified network does not have double edges between cells #
  
  # data.table variables
  unified_cells <- type_int <- N <- NULL
  
  spatial_network_annot <- GiottoUtils:::dt_sort_combine_two_columns(
    spatial_network_annot, "to", "from", "unified_cells"
  )
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
  
  # combine original and simulated network
  table_sim_results <- sample_dt[, .N, by = c(
    "unified_int", "type_int", "round"
  )]
  
  ## create complete simulations
  ## add 0 if no single interaction was found
  unique_ints <- unique(table_sim_results[, .(unified_int, type_int)])
  
  # data.table with 0's for all interactions
  minimum_simulations <- unique_ints[rep(
    seq_len(nrow(unique_ints)), number_of_simulations
  ), ]
  minimum_simulations[, round := rep(
    paste0("sim", seq_len(number_of_simulations)),
    each = nrow(unique_ints)
  )]
  minimum_simulations[, N := 0]
  
  table_sim_minimum_results <- rbind(table_sim_results, minimum_simulations)
  table_sim_minimum_results[, V1 := sum(N), by = c(
    "unified_int", "type_int", "round"
  )]
  table_sim_minimum_results <- unique(
    table_sim_minimum_results[, .(unified_int, type_int, round, V1)]
  )
  table_sim_results <- table_sim_minimum_results
  
  
  # data.table variables
  orig <- unified_int <- V1 <- original <- enrichm <- simulations <- NULL
  
  table_sim_results[, orig := "simulations"]
  spatial_network_annot[, round := "original"]
  
  table_orig_results <- spatial_network_annot[, .N, by = c(
    "unified_int", "type_int", "round"
  )]
  table_orig_results[, orig := "original"]
  data.table::setnames(table_orig_results, old = "N", new = "V1")
  
  table_results <- rbind(table_orig_results, table_sim_results)
  
  
  
  # add missing combinations from original or simulations
  # probably not needed anymore
  all_simulation_ints <- as.character(unique(table_results[
    orig == "simulations"
  ]$unified_int))
  all_original_ints <- as.character(unique(table_results[
    orig == "original"
  ]$unified_int))
  missing_in_original <- all_simulation_ints[
    !all_simulation_ints %in% all_original_ints
  ]
  missing_in_simulations <- all_original_ints[
    !all_original_ints %in% all_simulation_ints
  ]
  create_missing_for_original <- table_results[
    unified_int %in% missing_in_original
  ]
  create_missing_for_original <- unique(create_missing_for_original[
    , c("orig", "V1") := list("original", 0)
  ])
  create_missing_for_simulations <- table_results[
    unified_int %in% missing_in_simulations
  ]
  create_missing_for_simulations <- unique(
    create_missing_for_simulations[, c("orig", "V1") := list(
      "simulations", 0
    )]
  )
  # stop **********************************************************************************************
  
  table_results <- do.call(
    "rbind",
    list(
      table_results, create_missing_for_original,
      create_missing_for_simulations
    )
  )
  return(table_results)
}

cellProx.calcP <- function(table_results, number_of_simulations = 1000, 
                           adjust_method = c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel", "BY"), 
                           set_seed = TRUE, seed_number = 1234)
{
  sel_adjust_method <- match.arg(adjust_method, choices = c(
    "none", "fdr", "bonferroni", "BH",
    "holm", "hochberg", "hommel",
    "BY"
  ))
  ## p-values
  combo_list <- rep(NA, length = length(unique(table_results$unified_int)))
  p_high <- rep(NA, length = length(unique(table_results$unified_int)))
  p_low <- rep(NA, length = length(unique(table_results$unified_int)))
  
  for (int_combo in seq_along(unique(table_results$unified_int))) {
    this_combo <- as.character(unique(table_results$unified_int)[int_combo])
    
    sub <- table_results[unified_int == this_combo]
    
    orig_value <- sub[orig == "original"]$V1
    sim_values <- sub[orig == "simulations"]$V1
    
    length_simulations <- length(sim_values)
    if (length_simulations != number_of_simulations) {
      additional_length_needed <- number_of_simulations -
        length_simulations
      sim_values <- c(sim_values, rep(0, additional_length_needed))
    }
    
    p_orig_higher <- 1 - (sum((orig_value + 1) > (sim_values + 1)) /
                            number_of_simulations)
    p_orig_lower <- 1 - (sum((orig_value + 1) < (sim_values + 1)) /
                           number_of_simulations)
    
    combo_list[[int_combo]] <- this_combo
    p_high[[int_combo]] <- p_orig_higher
    p_low[[int_combo]] <- p_orig_lower
  }
  res_pvalue_DT <- data.table::data.table(
    unified_int = as.vector(combo_list),
    p_higher_orig = p_high,
    p_lower_orig = p_low
  )
  
  
  # depletion or enrichment in barplot format
  table_mean_results <- table_results[, .(mean(V1)), by = c(
    "orig", "unified_int", "type_int"
  )]
  table_mean_results_dc <- data.table::dcast.data.table(
    data = table_mean_results, formula = type_int + unified_int ~ orig,
    value.var = "V1"
  )
  table_mean_results_dc[, original := ifelse(is.na(original), 0, original)]
  table_mean_results_dc[, enrichm := log2((original + 1) / (simulations + 1))]
  
  
  table_mean_results_dc <- merge(
    table_mean_results_dc, res_pvalue_DT,
    by = "unified_int"
  )
  data.table::setorder(table_mean_results_dc, enrichm)
  table_mean_results_dc[, unified_int := factor(unified_int, unified_int)]
  
  # adjust p-values for mht
  
  # data.table variables
  p.adj_higher <- p.adj_lower <- p_lower_orig <- p_higher_orig <-
    PI_value <- int_ranking <- NULL
  
  table_mean_results_dc[, p.adj_higher := stats::p.adjust(
    p_higher_orig,
    method = sel_adjust_method
  )]
  table_mean_results_dc[, p.adj_lower := stats::p.adjust(
    p_lower_orig,
    method = sel_adjust_method
  )]
  
  
  table_mean_results_dc[, PI_value := ifelse(p.adj_higher <= p.adj_lower,
                                             -log10(p.adj_higher + (1 / number_of_simulations)) * enrichm,
                                             -log10(p.adj_lower + (1 / number_of_simulations)) * enrichm
  )]
  data.table::setorder(table_mean_results_dc, PI_value)
  
  # order
  table_mean_results_dc <- table_mean_results_dc[order(-PI_value)]
  table_mean_results_dc[, int_ranking := seq_len(.N)]
  
  return(list(
    raw_sim_table = table_results,
    enrichm_res = table_mean_results_dc
  ))
}


cellProx.subgroup <- function(dt)
{
  wide_dt <- dcast(dt,                # Your input data.table
                   unified_int ~ round, # Formula: rows ~ columns
                   value.var = "V1",  # Column with values for the new wide cells
                   fun.aggregate = mean, # Function to aggregate values (when multiple rows match a cell)
                   na.rm = TRUE,      # Tell mean() to ignore NA values during aggregation
                   fill = 0)    
  setkey(wide_dt, unified_int) # Good practice
  long_dt <- melt(wide_dt,
                  id.vars = "unified_int",       # Column(s) to keep as identifier variables
                  variable.name = "round",       # Name for the new column holding the former column names (1, 2, 3)
                  value.name = "V1",             # Name for the new column holding the values from the melted columns
                  variable.factor = FALSE)       # Keep the 'round' column as character initially (avoids factor issues)
  long_dt$orig <- "NA"
  long_dt$orig[grep("original", long_dt$round)] <- "original"
  long_dt$orig[grep("^sim", long_dt$round)] <- "simulations"
  long_dt$type_int <- "NA"
  long_dt$tosplit <- long_dt$unified_int
  long_dt <- as.data.table(tidyr::separate(long_dt, col = tosplit, into = c("left", "right"), sep = "--"))
  long_dt$type_int <- ifelse( long_dt$left == long_dt$right, "homo", "hetero")
  long_dt$left <- long_dt$right <- NULL
  out <- cellProx.calcP(long_dt, number_of_simulations = 200 )  # needs to match what was used to calculate sims
  return(out)
}


defineNeighborhoods <- function( annot.list)
{
  list_ID <- names(annot.list)
  res <- list()
  copymatrix <- list()
  for (i in list_ID)
  { 
    #                    from                   to sdimx_begin sdimy_begin sdimx_end sdimy_end distance     weight
    #                  <char>               <char>       <num>       <num>     <num>     <num>    <num>      <num>
    # 1: ctrl.1-ctrl.1_100000  ctrl.1-ctrl.1_99994      9482.0      9320.0    9482.5    9305.0 15.00833 0.06246747
    # 2: ctrl.1-ctrl.1_100000 ctrl.1-ctrl.1_100002      9482.0      9320.0    9500.0    9323.5 18.33712 0.05171401
    # 3: ctrl.1-ctrl.1_100000 ctrl.1-ctrl.1_100004      9482.0      9320.0    9462.5    9333.0 23.43608 0.04092309
    # 
    #                                      to_cell_type                           from_cell_type type_int
    #                                      <char>                                   <char>   <char>
    # 1:                    TCF1_MIXED_TOX_HI_CD8                       TCF1_HI_TOX_HI_CD8   hetero
    # 2:                           CD4_CD20_MIXED                       TCF1_HI_TOX_HI_CD8   hetero
    # 3:                             CD8_PC_MIXED                       TCF1_HI_TOX_HI_CD8   hetero
    # 
    #                                                    from_to                                              unified_int
    #                                                     <char>                                                   <char>
    # 1:                TCF1_HI_TOX_HI_CD8-TCF1_MIXED_TOX_HI_CD8                TCF1_HI_TOX_HI_CD8--TCF1_MIXED_TOX_HI_CD8
    # 2:                       TCF1_HI_TOX_HI_CD8-CD4_CD20_MIXED                       CD4_CD20_MIXED--TCF1_HI_TOX_HI_CD8
    # 3:                         TCF1_HI_TOX_HI_CD8-CD8_PC_MIXED                         CD8_PC_MIXED--TCF1_HI_TOX_HI_CD8
    
    print(paste0("pass 1: ", i))  # because the pairings are NOT DIRECTIONAL, I need to swap from/to and rbind it in
    a <- annot.list[[i]][, c(1:2,9:10)]  # simplify to just the necessary columns 
    b <- copy(a); names(b) <- c("to","from","from_cell_type", "to_cell_type")  # swap order of column names 
    a <- rbind(a, b)
    neighborhood.matrix <- a[, .N, by = .(from, to_cell_type, from_cell_type)] %>% 
      dcast(from + from_cell_type ~ to_cell_type, value.var= 'N', fill=0)  # count the cell types and make a wide data frame 
    neighborhood.matrix[, c('ARTIFACT', 'CD4.outlier', 'CD8.outlier'):= NULL]
    copymatrix[[i]] <- neighborhood.matrix
    
    #                    from                           from_cell_type ACTIVATED_CD4 ACTIVATED_CD8 ACTIVATED_MEMORY_CD4 
    #                  <char>                                   <char>         <int>         <int>                <int>              
    # 1: ctrl.1-ctrl.1_100000                       TCF1_HI_TOX_HI_CD8             0             0                    0                  
    # 2: ctrl.1-ctrl.1_100001                       TCF1_HI_TOX_HI_CD8             0             0                    0                  
    # 3: ctrl.1-ctrl.1_100002                           CD4_CD20_MIXED             0             0                    0  
    # 
    #    ACTIVATED_TREG ARTIFACT BCELL_MYELOID_MIXED BCELL_MYELOID_NONIMMUNE_MIXED CD4.outlier CD4_CD20_MIXED CD4_CD8_CD20_MIXED
    #             <int>    <int>               <int>                         <int>       <int>          <int>              <int>
    # 1:              0        0                   0                             1           0              1                  2
    # 2:              0        0                   0                             0           0              1                  0
    # 3:              0        0                   0                             1           0              0                  2
  }
  for (i in list_ID)
  { 
    print(paste0("pass 2: ", i))
    working.dt <- copy(copymatrix[[i]])
    id_cols <- c("from", "from_cell_type")
    count_cols <- setdiff(names(working.dt), id_cols)
    
    # 2. Calculate the total for each row and add it as a new column
    working.dt[, row_total := rowSums(.SD), .SDcols = count_cols]
    
    # 3. Divide each count column by the row_total to get the proportion
    #    We loop through the count columns and update them "by reference"
    working.dt[, (count_cols) := lapply(.SD, function(x) x / row_total), .SDcols = count_cols]
    res[[i]] <- working.dt
  } 
  return( list(prop.table=res, raw.table=copymatrix, list_ID=list_ID)) 
}



niche.contrast <- function( annot.table, indexCluster, common.title, grouping.exp, grouping.ctrl, condition.exp, condition.ctrl)
{
  ctrl <- showNiche.multiSample(annot.table[grouping.ctrl], indexCluster, title = paste0(common.title,".",condition.ctrl),
                                plots=F)
  expt <- showNiche.multiSample(annot.table[grouping.exp], indexCluster, title = paste0(common.title,".",condition.exp),
                                plots=F)
  ctrl$stage <- condition.ctrl;  expt$stage <- condition.exp
  res <- rbind(ctrl, expt)
  res$stage <- factor(res$stage, levels = c(condition.exp, condition.ctrl))
  res$subset <- factor(res$subset, levels = levels(expt$subset))
  stableColorsL2.fxn <- stableColorsL2
  
  # Create a lighter version of the color palette
  lighter_stableColorsL2.fxn <- sapply(stableColorsL2.fxn, lighten_hex, factor = 0.6) # Adjust factor as needed
  
  # Create a combined color mapping key
  res$color_key <- interaction(res$subset, res$stage, sep = "_")
  
  # Create the combined color palette
  combined_palette <- c(stableColorsL2.fxn, lighter_stableColorsL2.fxn[names(stableColorsL2.fxn)])
  names(combined_palette) <- c(paste0(names(stableColorsL2.fxn), "_", condition.exp), 
                               paste0(names(lighter_stableColorsL2.fxn), "_", condition.ctrl))
  
  annotateX <- max(res$proportion)-0.2*max(res$proportion)
  print(a <- ggplot(res, aes(y = subset, x = proportion, fill=color_key, label=list_ID, group=stage)) +
          geom_point(size=2, position=position_dodge(width=1), pch=21, stroke=0.2, color='grey50') +
          annotate('text',x= annotateX ,y=35, label=paste0(condition.ctrl," in light shade\n", 
                                                           condition.exp," in dark shade"), size=3)+
          stat_summary( fun.min = "mean", fun.max = "mean", geom = "errorbar", position = position_dodge(width=1),
                        width = 1, color = "black", linewidth = 0.5 ) +
          theme_bw() + ggtitle(common.title) +
          scale_fill_manual(name = 'subset', values= combined_palette,  guide = guide_legend(ncol = 2))  +
          theme(axis.text = element_text(size=6), axis.title.y = element_blank()) +
          theme(legend.position = 'none')
  )
  return(list(res=res, scatterplot = a, colorScheme = combined_palette))
}


















permutation.test.neighborhoods <- function(num.permutations = 10000, data.tbl, group.0, group.1)
{
  obs.diff <- mean(data.tbl[source %in% group.0, data]) - mean(data.tbl[source %in% group.1, data])  
  print(paste0("Observed difference: ", round(obs.diff, 2)))  
  permutation.diffs <- numeric(num.permutations)  # now for a ridiculous number of permutations 
  for (i in 1:num.permutations) 
  {
    shuffled_source <- sample(data.tbl$source)    # randomly shuffle the group labels.
    
    # Calculate the difference in means for this random permutation.
    # 1. mutate(): Replaces the original 'source' column with the 'shuffled_source' labels.
    # 2. group_by(): Groups the data by the new, random source labels.
    # 3. summarise(): Calculates the mean for each of the two groups.
    # 4. pull(): Extracts the mean values as a vector.
    # 5. diff(): Calculates the difference between the two means in the vector.
    permutation.diffs[i] <- data.tbl %>%
      mutate(source = shuffled_source) %>%
      group_by(source) %>%
      summarise(mean = mean(data)) %>%
      pull(mean) %>%
      diff()
  }
  hist(permutation.diffs, breaks=100)
  # The p-value is the proportion of permuted differences that were as extreme as (or more extreme than) our observed difference. 
  # The 'abs()' function performs a two-tailed test. The logical operator yields T/F and then mean calculates the proportion 
  pV <- mean(abs(permutation.diffs) >= abs(obs.diff))
  print(paste0("P value is ", pV, ". Lowest possible P value for ",num.permutations, " iterations is ", 1/num.permutations))
  return(pV)
}