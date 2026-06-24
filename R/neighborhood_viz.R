#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarise ungroup mutate arrange select desc
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
#' @return A ggplot object or a data.table summary tracking distribution properties.
#' @export
showNicheSingleSample <- function(annot.table,
                                  indexCluster,
                                  title = NULL,
                                  plots = TRUE,
                                  color_vector = NULL,
                                  plot_type = "lollipop") {

  res <- annot.table[!unified_int %like% "ARTIFACT" & !unified_int %like% "outlier"]

  res <- rbind(
    res[from_cell_type == indexCluster, .(cell_type = to_cell_type)],
    res[to_cell_type == indexCluster, .(cell_type = from_cell_type)]
  )

  freq_table <- table(res$cell_type)
  prop_table <- prop.table(freq_table)
  res_dt <- data.table::as.data.table(prop_table)
  colnames(res_dt) <- c("subset", "proportion")

  res_dt <- res_dt[order(-proportion)]
  res_dt[, subset := factor(subset, levels = subset)]

  if (!plots) return(res_dt)
  if (is.null(title)) title <- indexCluster

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
#' @param annot.list A list of annotated spatial network tables.
#' @param indexCluster Character string. The "center" cell type.
#' @param title Character string. Title for the plot.
#' @param plots Logical. Returns plot if TRUE.
#' @param color_vector Named character vector of colors.
#'
#' @return A ggplot object or summary reference table framework.
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
#' @return A ggraph visualization layout framework object.
#' @export
plotAverageNicheIgraph <- function(annot.list,
                                   indexCluster,
                                   title = NULL,
                                   color_vector = NULL,
                                   top_n = 8) {

  res <- showNicheMultiSample(annot.list, indexCluster, plots = FALSE)
  ave <- res %>%
    dplyr::group_by(subset) %>%
    dplyr::summarise(mean_proportion = mean(proportion, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(mean_proportion)) %>%
    head(top_n)

  center_node <- data.frame(subset = indexCluster, mean_proportion = 0)
  grid_data <- rbind(center_node, as.data.frame(ave))
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
    p := p + scale_fill_manual(values = color_vector, guide = "none")
  }

  return(p)
}

#' Contrast Niche Composition Between Conditions
#'
#' Visualizes the difference in niche composition for a specific cell type between two experimental groups.
#'
#' @param annot.list A list of annotated spatial network tables.
#' @param indexCluster Character string. The center cell type.
#' @param group.exp Character vector of sample IDs for the experimental group.
#' @param group.ctrl Character vector of sample IDs for the control group.
#' @param condition.exp Label for experimental group.
#' @param condition.ctrl Label for control group.
#' @param color_vector Named character vector of colors.
#'
#' @return A list containing the ggplot object (`scatterplot`) and the raw calculation table (`res`).
#' @export
nicheContrast <- function(annot.list,
                          indexCluster,
                          group.exp,
                          group.ctrl,
                          condition.exp = "Exp",
                          condition.ctrl = "Ctrl",
                          color_vector = NULL) {

  ctrl <- showNicheMultiSample(annot.list[group.ctrl], indexCluster, plots = FALSE)
  ctrl[, stage := condition.ctrl]

  expt <- showNicheMultiSample(annot.list[group.exp], indexCluster, plots = FALSE)
  expt[, stage := condition.exp]

  res <- rbind(ctrl, expt)
  res[, stage := factor(stage, levels = c(condition.exp, condition.ctrl))]

  subset_levels <- expt %>%
    dplyr::group_by(subset) %>%
    dplyr::summarise(m = mean(proportion)) %>%
    dplyr::arrange(dplyr::desc(m)) %>%
    dplyr::pull(subset)

  res[, subset := factor(subset, levels = subset_levels)]

  p <- ggplot(res, aes(y = subset, x = proportion, fill = subset, group = stage)) +
    geom_point(size = 2, position = position_dodge(width = 0.75), pch = 21, stroke = 0.2, color = "grey50") +
    stat_summary(fun.min = mean, fun.max = mean, geom = "errorbar",
                 position = position_dodge(width = 0.75), width = 0.5) +
    theme_bw() +
    ggtitle(paste(indexCluster, ":", condition.exp, "vs", condition.ctrl)) +
    theme(legend.position = "none")

  if (!is.null(color_vector)) {
    p <- p + scale_fill_manual(values = color_vector)
  }

  return(list(scatterplot = p, res = res))
}

#' Plot Interaction Effect Sizes
#'
#' Visualizes calculated cell-cell affinity effect size metrics as a dot/bar plot.
#'
#' @param df Data.table output from `calculateInteractionEffectSizes`.
#' @param condition.exp Label for experimental condition.
#' @param condition.ctrl Label for control condition.
#'
#' @return A ggplot object.
#' @export
plotInteractionEffectSizes <- function(df, condition.exp, condition.ctrl) {
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

#' Plot Neighbor-Dependent Feature Expression Histogram
#'
#' Compares feature expression in cells with neighbors of type X vs those without.
#'
#' @param gobject A Giotto object.
#' @param neighborhood_prop_table Data.table from `defineNeighborhoods` (one sample).
#' @param target_cluster The cell type of the "center" cell.
#' @param neighbor_cluster The neighbor cell type to stratify by.
#' @param feature Gene or protein name to plot.
#' @param title Plot title.
#'
#' @return A ggplot rendering framework.
#' @export
plotNeighborHistogram <- function(gobject,
                                  neighborhood_prop_table,
                                  target_cluster,
                                  neighbor_cluster,
                                  feature,
                                  title = NULL) {

  target_cells <- neighborhood_prop_table[from_cell_type == target_cluster]
  cells_with <- target_cells[get(neighbor_cluster) > 0, from]
  cells_without <- target_cells[get(neighbor_cluster) == 0, from]

  if (length(cells_with) == 0 || length(cells_without) == 0) {
    warning("One group has 0 cells. Cannot plot.")
    return(NULL)
  }

  expr_with <- Giotto::getExpression(gobject, values = "normalized")[feature, cells_with]
  expr_without := Giotto::getExpression(gobject, values = "normalized")[feature, cells_without]

  dt <- data.table(
    value = c(expr_with, expr_without),
    group = factor(c(rep("Neighbors", length(expr_with)), rep("No Neighbors", length(expr_without))),
                   levels = c("No Neighbors", "Neighbors"))
  )

  p <- ggplot(dt, aes(x = value, fill = group)) +
    geom_histogram(binwidth = 0.25, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = c("grey50", "firebrick")) +
    theme_classic() +
    ggtitle(ifelse(is.null(title), paste(feature, "in", target_cluster), title)) +
    xlab(paste(feature, "Expression")) +
    theme(legend.position = "top")

  return(p)
}

#' Plot Global Effect Size Metrics as Violin Layouts
#'
#' Evaluates and renders pooled delta shift distribution profiles across cell lineages.
#'
#' @param importIntxn.obj Comprehensive output object compiled via \code{findImportantInteractions}.
#' @param colNames String vector designating cell subset categories to track.
#' @param title Custom text heading to display over plot window layout.
#'
#' @return A ggplot visualization plot framework.
#' @export
plot.globalChanges.Effectsize <- function(importIntxn.obj, colNames, title = "Stage 3 vs Stage 1") {
  raw_data <- as.data.table(importIntxn.obj$fullggplot$data)
  res_list <- lapply(X = colNames, FUN = function(x) { raw_data[unified_int %like% x, .(effect.size)] })
  names(res_list) <- colNames

  b <- purrr::list_rbind(res_list, names_to = 'subset')
  names(b)[2] <- "Effect.size"
  b <- b[!is.nan(Effect.size) & !subset %like% 'ARTIFACT' & !subset %like% 'outlier' & !subset %like% 'Tumor']

  order_dt <- b %>% dplyr::group_by(subset) %>% dplyr::summarise(mean_val = mean(Effect.size, na.rm = TRUE)) %>% dplyr::arrange(mean_val)
  b$subset <- factor(b$subset, levels = order_dt$subset)

  p <- ggplot(b, aes(y = subset, x = Effect.size)) +
    geom_vline(xintercept = 0, color = 'grey50') +
    geom_violin(fill = "grey90") +
    theme_bw() +
    theme(axis.title.y = element_blank(), axis.text = element_text(color = 'black')) +
    ggtitle(title)
  return(p)
}

#' Generate Comprehensive Rank Order Effect Size Scatter Plots
#'
#' Constructs clear distribution layouts tracking microenvironment perturbation shifts.
#'
#' @param imp.intxn.object Tracking framework returned via \code{findImportantInteractions}.
#' @param ef.threshold Cutoff constraint isolating relevant target metrics rows.
#' @param leaveOutClusters String keys identifying noise structures to filter out from rendering view.
#' @param labelSubset String indicators isolating specific subsets to attach text annotations onto.
#'
#' @return A structured rank distribution plot framework.
#' @export
plot.interaction.rankOrder <- function(imp.intxn.object, ef.threshold = 2, leaveOutClusters = NULL, labelSubset = NULL) {
  df <- data.table::copy(imp.intxn.object[[1]])
  if(!is.null(leaveOutClusters)) {
    for (cl in leaveOutClusters) { df <- df[!unified_int %like% cl] }
  }
  df <- df[abs(effect.size) > ef.threshold]
  df[, rankOrder := .I]
  df[, labelSubset := ""]

  max_rank <- max(df$rankOrder, na.rm = TRUE)
  min_rank <- min(df$rankOrder, na.rm = TRUE)
  annotateX.ctrl <- max_rank * 0.95
  annotateX.exp <- min_rank + 0.2 * min_rank
  annotateY <- min(df$effect.size, na.rm = TRUE) + 0.2 * min(df$effect.size, na.rm = TRUE)

  if (!is.null(labelSubset)) {
    indices_to_label <- rep(FALSE, nrow(df))
    for (lb in labelSubset) { indices_to_label <- indices_to_label | (df$unified_int %like% lb) }
    df[indices_to_label, labelSubset := as.character(unified_int)]
    df[, labelSubset := gsub("--", " \u2022 ", labelSubset)]
    p <- ggplot(df, aes(x = rankOrder, y = effect.size, label = labelSubset))
  } else {
    p <- ggplot(df, aes(x = rankOrder, y = effect.size, label = unified_int))
  }

  df[, target_nudge_x := ifelse(effect.size >= 0, max_rank * 0.25, -max_rank * 0.25)]
  df[, target_nudge_y := ifelse(effect.size >= 0, max(effect.size) * 0.15, min(effect.size) * 0.15)]

  nudge_x_vector <- if(!is.null(labelSubset)) ifelse(df$labelSubset != "", df$target_nudge_x, 0) else df$target_nudge_x
  nudge_y_vector := if(!is.null(labelSubset)) ifelse(df$labelSubset != "", df$target_nudge_y, 0) else df$target_nudge_y

  p <- p + geom_point(pch = 21, color = 'black', fill = 'gray50', alpha = 0.3, size = 5) +
    ggrepel::geom_text_repel(
      max.overlaps = 5, size = 3,
      nudge_x = nudge_x_vector, nudge_y = nudge_y_vector,
      direction = "both", segment.color = "grey30",
      box.padding = 0.8, point.padding = 0.4, force = 10, min.segment.length = 0
    ) +
    geom_hline(yintercept = 0, color = 'black') +
    ggtitle(paste0("Effect size: ", imp.intxn.object$parameters$condition.exp, " vs ", imp.intxn.object$parameters$condition.ctrl)) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 8, color = 'black'), plot.title = element_text(size = 10),
          axis.title.x = element_text(size = 10, color = 'black'), axis.text.x = element_text(size = 10, color = 'black')) +
    theme(legend.position = 'none') +
    annotate("text", x = annotateX.exp, y = annotateY, size = 5, label = imp.intxn.object$parameters$condition.exp) +
    annotate("text", x = annotateX.ctrl, y = annotateY, size = 3, label = imp.intxn.object$parameters$condition.ctrl) +
    ylab("Effect size") + coord_cartesian(clip = 'off')

  return(p)
}
