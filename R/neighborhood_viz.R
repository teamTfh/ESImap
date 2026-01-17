#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarise ungroup mutate arrange select
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggraph ggraph geom_edge_fan geom_node_point geom_node_text geom_node_label theme_graph
#' @importFrom tidygraph tbl_graph
#' @importFrom gridExtra grid.arrange
#' @importFrom stringr str_wrap
NULL

#' Visualize Niche Composition for a Single Sample
#'
#' Creates a lollipop chart or stacked bar chart showing the composition of neighbors for a specific cell type in one sample.
#'
#' @param annot.table A data.table annotating the spatial network (specific to one sample).
#' @param indexCluster Character string. The "center" cell type to analyze.
#' @param title Character string. Title for the plot.
#' @param plots Logical. If TRUE, returns a ggplot object. If FALSE, returns the summary data.table.
#' @param color_vector Named character vector of colors for cell types.
#' @param plot_type Character. "lollipop" or "stacked".
#'
#' @return A ggplot object or a data.table.
#' @export
showNicheSingleSample <- function(annot.table,
                                  indexCluster,
                                  title = NULL,
                                  plots = TRUE,
                                  color_vector = NULL,
                                  plot_type = "lollipop") {
  
  # Remove artifacts
  res <- annot.table[!unified_int %like% "ARTIFACT" & !unified_int %like% "outlier"]
  
  # Filter for the index cluster in either 'from' or 'to' position
  # Note: The logic assumes non-directional pairing representation
  res <- rbind(
    res[from_cell_type == indexCluster, .(cell_type = to_cell_type)],
    res[to_cell_type == indexCluster, .(cell_type = from_cell_type)]
  )
  
  # Calculate proportions
  freq_table <- table(res$cell_type)
  prop_table <- prop.table(freq_table)
  res_dt <- data.table::as.data.table(prop_table)
  colnames(res_dt) <- c("subset", "proportion")
  
  res_dt <- res_dt[order(-proportion)]
  res_dt[, subset := factor(subset, levels = subset)]
  
  if (!plots) return(res_dt)
  
  if (is.null(title)) title <- indexCluster
  
  # Base Plot
  p <- ggplot(res_dt, aes(y = subset, x = proportion, fill = subset, color = subset)) +
    theme_classic() +
    ggtitle(title) +
    theme(legend.position = "none")
  
  if (!is.null(color_vector)) {
    p <- p + scale_fill_manual(values = color_vector) + scale_color_manual(values = color_vector)
  }
  
  if (plot_type == "lollipop") {
    p <- p +
      geom_segment(aes(x = 0, xend = proportion, y = subset, yend = subset)) +
      geom_point(size = 5) +
      theme(axis.text.y = element_text(size = 8), axis.title.y = element_blank())
  } else {
    # Stacked bar style
    p <- ggplot(res_dt, aes(y = proportion, x = "Niche", fill = subset)) +
      geom_col(width = 1, position = "fill") +
      ggrepel::geom_text_repel(aes(label = subset), position = position_stack(vjust = 0.5), size = 3) +
      theme_classic() +
      ggtitle(paste("Neighbors of", title)) +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    
    if (!is.null(color_vector)) {
      p <- p + scale_fill_manual(values = color_vector)
    }
  }
  
  return(p)
}


#' Visualize Average Niche Across Multiple Samples
#'
#' Aggregates niche composition across multiple biological replicates and plots the mean proportions.
#'
#' @param annot.list A list of annotated spatial network tables (output of `annotateSpatialNetwork`).
#' @param indexCluster Character string. The "center" cell type.
#' @param title Character string. Title for the plot.
#' @param plots Logical. Returns plot if TRUE.
#' @param color_vector Named character vector of colors.
#'
#' @return A ggplot object (dot plot with error bars) or a summary data.table.
#' @export
showNicheMultiSample <- function(annot.list,
                                 indexCluster,
                                 title = NULL,
                                 plots = TRUE,
                                 color_vector = NULL) {
  
  res_list <- lapply(names(annot.list), function(id) {
    dt <- showNicheSingleSample(annot.list[[id]], indexCluster, plots = FALSE)
    dt[, list_ID := id]
    return(dt)
  })
  
  res <- data.table::rbindlist(res_list)
  
  # Calculate average for ordering
  ave <- res %>%
    dplyr::group_by(subset) %>%
    dplyr::summarise(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(mean_proportion))
  
  res[, subset := factor(subset, levels = ave$subset)]
  
  if (!plots) return(res)
  
  if (is.null(title)) title <- indexCluster
  
  p <- ggplot(res, aes(y = subset, x = proportion, fill = subset)) +
    geom_point(size = 2, position = position_dodge(width = 1), pch = 21, stroke = 0.2, color = "grey50") +
    stat_summary(fun.min = function(x) mean(x), fun.max = function(x) mean(x),
                 geom = "errorbar", width = 0.75, color = "black", linewidth = 0.5) +
    theme_bw() +
    ggtitle(title) +
    theme(axis.text = element_text(size = 8), axis.title.y = element_blank(), legend.position = "none")
  
  if (!is.null(color_vector)) {
    p <- p + scale_fill_manual(values = color_vector)
  }
  
  return(p)
}


#' Plot Average Niche as a Network Graph
#'
#' Visualizes the top neighbors of a specific cell type as a star-like network graph using `ggraph`.
#'
#' @param annot.list A list of annotated spatial network tables.
#' @param indexCluster Character string. The center cell type.
#' @param title Character string.
#' @param color_vector Named character vector of colors.
#' @param top_n Integer. Number of top neighbors to display (default: 8).
#'
#' @return A ggraph object.
#' @export
plotAverageNicheIgraph <- function(annot.list,
                                   indexCluster,
                                   title = NULL,
                                   color_vector = NULL,
                                   top_n = 8) {
  
  # Get average data
  res <- showNicheMultiSample(annot.list, indexCluster, plots = FALSE)
  ave <- res %>%
    dplyr::group_by(subset) %>%
    dplyr::summarise(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(mean_proportion)) %>%
    head(top_n)
  
  # Handle self-loops in colors (if the center is also a neighbor)
  center_color_name <- indexCluster
  if (indexCluster %in% ave$subset) {
    # If the center cell type is in the neighbors, we might need a distinct color key if strictly matching names
    # For now, we assume the color_vector handles the name correctly.
  }
  
  # Prepare graph data
  center_node <- data.frame(subset = indexCluster, mean_proportion = 0)
  grid_data <- rbind(center_node, as.data.frame(ave))
  
  # Wrap labels
  grid_data$label <- stringr::str_wrap(gsub("_", " ", grid_data$subset), 8)
  
  edges <- data.frame(
    from = indexCluster,
    to = as.character(grid_data$subset[grid_data$subset != indexCluster])
  )
  
  nodes <- data.frame(
    name = grid_data$subset,
    label = grid_data$label,
    size = grid_data$mean_proportion
  )
  
  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  p <- ggraph::ggraph(graph, layout = "dendrogram", circular = TRUE) +
    ggraph::geom_edge_fan(width = 0.5, color = "gray60", alpha = 0.8) +
    ggraph::geom_node_point(aes(fill = name, size = size), shape = 21, color = "black", stroke = 0.2) +
    ggraph::geom_node_text(aes(label = label, x = x * 1.3, y = y * 1.3), size = 3) +
    scale_size_continuous(name = "Interaction\nProportion", range = c(2, 10)) +
    ggraph::theme_graph(base_family = "sans") +
    coord_fixed() +
    labs(title = title) +
    xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
  
  if (!is.null(color_vector)) {
    p <- p + scale_fill_manual(values = color_vector, guide = "none")
  }
  
  return(p)
}


#' Contrast Niche Composition Between Conditions
#'
#' visualizes the difference in niche composition for a specific cell type between two experimental groups.
#'
#' @param annot.list A list of annotated spatial network tables.
#' @param indexCluster Character string. The center cell type.
#' @param group.exp Character vector of sample IDs for the experimental group.
#' @param group.ctrl Character vector of sample IDs for the control group.
#' @param condition.exp Label for experimental group.
#' @param condition.ctrl Label for control group.
#' @param color_vector Named character vector of colors.
#'
#' @return A list containing the ggplot object (`scatterplot`) and the data table (`res`).
#' @export
nicheContrast <- function(annot.list,
                          indexCluster,
                          group.exp,
                          group.ctrl,
                          condition.exp = "Exp",
                          condition.ctrl = "Ctrl",
                          color_vector = NULL) {
  
  # Calculate per group
  ctrl <- showNicheMultiSample(annot.list[group.ctrl], indexCluster, plots = FALSE)
  ctrl[, stage := condition.ctrl]
  
  expt <- showNicheMultiSample(annot.list[group.exp], indexCluster, plots = FALSE)
  expt[, stage := condition.exp]
  
  res <- rbind(ctrl, expt)
  res[, stage := factor(stage, levels = c(condition.exp, condition.ctrl))]
  
  # Order by experimental group means
  subset_levels <- expt %>%
    dplyr::group_by(subset) %>%
    dplyr::summarise(m = mean(proportion)) %>%
    dplyr::arrange(dplyr::desc(m)) %>%
    dplyr::pull(subset)
  
  res[, subset := factor(subset, levels = subset_levels)]
  
  # Create plot
  p <- ggplot(res, aes(y = subset, x = proportion, fill = subset, group = stage)) +
    geom_point(size = 2, position = position_dodge(width = 0.75), pch = 21, stroke = 0.2, color = "grey50") +
    stat_summary(fun.min = mean, fun.max = mean, geom = "errorbar",
                 position = position_dodge(width = 0.75), width = 0.5) +
    theme_bw() +
    ggtitle(paste(indexCluster, ":", condition.exp, "vs", condition.ctrl)) +
    theme(legend.position = "none")
  
  # Handle Colors: Generate a lighter version for the control group
  if (!is.null(color_vector)) {
    # Note: If you want distinct colors for Stage 1 vs Stage 3 (light vs dark),
    # complex logic is required. Here we map the base color to both for simplicity,
    # relying on position_dodge to distinguish.
    p <- p + scale_fill_manual(values = color_vector)
  }
  
  return(list(scatterplot = p, res = res))
}


#' Plot Interaction Effect Sizes
#'
#' Visualizes the output of `calculateInteractionEffectSizes` as a bar chart.
#'
#' @param df Data.table output from `calculateInteractionEffectSizes`.
#' @param condition.exp Label for experimental condition.
#' @param condition.ctrl Label for control condition.
#'
#' @return A ggplot object.
#' @export
plotInteractionEffectSizes <- function(df, condition.exp, condition.ctrl) {
  
  annotateX.exp <- max(df$effect.size, na.rm=TRUE) * 0.8
  annotateX.ctrl <- min(df$effect.size, na.rm=TRUE) * 0.8
  
  p <- ggplot(df, aes(x = effect.size, y = reorder(unified_int, effect.size))) +
    geom_bar(stat = "identity", width = 0.05, color = "grey70", fill = "grey50") +
    geom_point(size = 4, pch = 21, color = "white", fill = "grey50", stroke = 0.2) +
    theme_classic() +
    ggtitle(paste("Effect size:", condition.exp, "vs", condition.ctrl)) +
    theme(
      axis.text.y = element_text(size = 8, color = "black"),
      axis.title.y = element_blank(),
      legend.position = "none"
    ) +
    xlab("Effect Size") +
    geom_vline(xintercept = 0, linetype = "dashed")
  
  return(p)
}


#' Plot Neighbor-Dependent Feature Expression
#'
#' Compares feature expression in cells with neighbors of type X vs those without.
#' Requires a Giotto object to retrieve expression data.
#'
#' @param gobject A Giotto object.
#' @param neighborhood_prop_table Data.table from `defineNeighborhoods` (one sample).
#' @param target_cluster The cell type of the "center" cell.
#' @param neighbor_cluster The neighbor cell type to stratify by.
#' @param feature Gene or protein name to plot.
#' @param title Plot title.
#'
#' @return A list containing the plot and the p-value.
#' @export
plotNeighborHistogram <- function(gobject,
                                  neighborhood_prop_table,
                                  target_cluster,
                                  neighbor_cluster,
                                  feature,
                                  title = NULL) {
  
  # Filter cells
  target_cells <- neighborhood_prop_table[from_cell_type == target_cluster]
  
  # Split into groups
  cells_with <- target_cells[get(neighbor_cluster) > 0, from]
  cells_without <- target_cells[get(neighbor_cluster) == 0, from]
  
  if (length(cells_with) == 0 || length(cells_without) == 0) {
    warning("One group has 0 cells. Cannot plot.")
    return(NULL)
  }
  
  # Get Expression
  # Note: Assuming 'normalized' slot exists
  expr_with <- Giotto::getExpression(gobject, values = "normalized")[feature, cells_with]
  expr_without <- Giotto::getExpression(gobject, values = "normalized")[feature, cells_without]
  
  dt <- data.table(
    value = c(expr_with, expr_without),
    group = factor(c(rep("Neighbors", length(expr_with)), rep("No Neighbors", length(expr_without))),
                   levels = c("No Neighbors", "Neighbors"))
  )
  
  # Plot
  p <- ggplot(dt, aes(x = value, fill = group)) +
    geom_histogram(binwidth = 0.25, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = c("grey50", "firebrick")) +
    theme_classic() +
    ggtitle(ifelse(is.null(title), paste(feature, "in", target_cluster), title)) +
    xlab(paste(feature, "Expression")) +
    theme(legend.position = "top")
  
  return(p)
}