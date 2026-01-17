
rgb_to_hex <- function(rgb) {
  sprintf("#%02X%02X%02X", rgb[1], rgb[2], rgb[3])
}

lighten_hex <- function(hex, factor = 0.3) { # You can adjust the factor
  rgb <- col2rgb(hex)
  new_rgb <- round(rgb + (255 - rgb) * factor)
  new_rgb <- pmax(0, pmin(255, new_rgb))
  return(rgb_to_hex(new_rgb))
}


removeArtifacts <- function(dt)
{
  dt <- dt[!unified_int %like% "ARTIFACT"]
  dt <- dt[!unified_int %like% "outlier"]
  return(dt)
}



showNiche.singleSample <- function(annot.table, indexCluster, title, plots=TRUE)
{
  # example annot.table[['tumor.5']]
  #                        from                     to sdimx_begin sdimy_begin sdimx_end sdimy_end distance     weight
  #                    <char>                 <char>       <num>       <num>     <num>     <num>    <num>      <num>
  # 1: tumor.5-tumor.5_100000  tumor.5-tumor.5_99987     39429.5     14436.0   39446.5   14427.0 19.23538 0.04941838
  # 2: tumor.5-tumor.5_100000  tumor.5-tumor.5_99971     39429.5     14436.0   39428.5   14416.5 19.52562 0.04871959
  # 3: tumor.5-tumor.5_100000  tumor.5-tumor.5_99999     39429.5     14436.0   39408.0   14439.0 21.70829 0.04403677
  #                    to_cell_type              from_cell_type type_int                                                 from_to
  #                         <char>                      <char>   <char>                                                  <char>
  # 1:                  MEMORY_CD4                  MEMORY_CD8   hetero                                   MEMORY_CD8-MEMORY_CD4
  # 2:                  MEMORY_CD4                  MEMORY_CD8   hetero                                   MEMORY_CD8-MEMORY_CD4
  # 3:                ICOS_LO_TREG                  MEMORY_CD8   hetero                                 MEMORY_CD8-ICOS_LO_TREG
  #                                                       unified_int
  #                                                      <char>
  # 1:                                   MEMORY_CD4--MEMORY_CD8
  # 2:                                   MEMORY_CD4--MEMORY_CD8
  # 3:                                 ICOS_LO_TREG--MEMORY_CD8
  # 
  
  res <- removeArtifacts(annot.table)
  res <- rbind( res[from_cell_type == indexCluster, to_cell_type], 
                res[to_cell_type == indexCluster, from_cell_type] )  
  # because the from/to is not really directional, and we'd be missing some of the parings if we assume directionality
  res <- res %>% table %>% prop.table  # calculate proportions 
  res <- as.data.table(res); colnames(res) <- c('subset','proportion')
  res <- res[order(res$proportion, decreasing=T)]
  stableColorsL2.fxn <- stableColorsL2
  
  res$subset <- factor(res$subset, levels = res$subset)   
  
  if(plots == T)
  {  
    # lollipop chart 
    print( 
      a <- ggplot(data = res, aes(y = subset, x = proportion, fill=subset, color = subset)) + 
        geom_bar(stat='identity', width=0.01) + 
        geom_point(size=5) + 
        theme_classic() + ggtitle(title) + 
        scale_fill_manual(values = stableColorsL2.fxn) +
        scale_color_manual(values = stableColorsL2.fxn) +
        theme(axis.text.y = element_text( size=6), axis.title.y = element_blank(), plot.title=element_text(size=8)) + 
        theme(legend.position = 'none')  
    )
    
    # stacked bar chart 
    print(
      ggplot(data = res, aes(y = proportion, x = 'blank',  fill = subset)) + 
        geom_col(width=1, position="fill") + 
        ggrepel::geom_text_repel( aes(label = subset ),
                                  position = position_stack(vjust = 0.5),           # specify stacked position to follow bar layer's default position
                                  size = 2,  color='black', direction = "y", xlim = c(1.5, NA), hjust = 0,
                                  segment.linetype = "dotted", box.padding = .4 ) +
        theme_classic() + ggtitle(paste("Neighbors of ", title)) +  
        theme(axis.text.x = element_blank()) +
        theme(axis.title.x = element_blank()) +
        # theme(legend.text = element_text(size=3), legend.title = element_text(size=3)) + 
        theme(legend.position = 'none') +
        # theme( plot.margin = margin(0,200,0,0,  unit="pt")) + 
        scale_fill_manual(values = stableColorsL2.fxn) + scale_color_manual(values = stableColorsL2.fxn) 
    )
  }
  if(plots==T)   return(a)
  else return(res)
}



showNiche.multiSample <- function(annot.table, indexCluster, title,  plots=T, parallel=F)
{
  library(future.apply)
  res <- list()
  IDlist <- names(annot.table)
  if(parallel == T)
  {  
    res_list <- future_lapply(IDlist, function(i) {
      single_res <- showNiche.singleSample( annot.table = annot.table[[i]], indexCluster = indexCluster, title = title, plots = F )
      single_res$list_ID <- i    # name by list_ID
      return(single_res)
    }, future.seed = TRUE) # future.seed = TRUE is crucial for reproducible results!
  }  
  
  if (parallel == F)
  {  
    for (i in IDlist)
    {
      res[[i]] <- showNiche.singleSample(annot.table = annot.table[[i]], indexCluster = indexCluster,
                                         title=title, plots=F )  # return per-sample proportions
      res[[i]]$list_ID <- i  # name by list_ID
    }
  }
  res <- rbindlist(res)  # merge into one data table 
  ave <-  res %>% group_by(subset) %>% summarise(mean_proportion = mean(proportion, na.rm = TRUE))  # calculate mean over grouping
  ave <- ave[order(ave$mean_proportion, decreasing = T),]  # order by decreasing value 
  res$subset <- factor(res$subset, levels = ave$subset)  # set the res$subset factors based on the decreasing average
  stableColorsL2.fxn <- stableColorsL2
  
  if(plots == T)
  { 
    # scatter plot of biological replicates
    print( 
      ggplot(res, aes(y = subset, x = proportion, fill=subset, label=list_ID)) +
        geom_point(size=2, position=position_dodge(width=1), pch=21, stroke=0.2, color='grey50') +
        stat_summary( fun.min = "mean", fun.max = "mean", geom = "errorbar",
                      width = 0.75, color = "black", linewidth = 0.5 ) + # geom_path(aes(group = list_ID), alpha=0.2, color='grey80') + 
        theme_bw() + ggtitle(title) + 
        # ggrepel::geom_label_repel(size=2, box.padding = unit(0.1,"lines"), alpha=0.2) + 
        scale_fill_manual(values= stableColorsL2.fxn)  + 
        theme(axis.text = element_text(size=6), axis.title.y = element_blank()) + 
        theme(legend.position = 'none')
    )
  }
  return(res)
}  



plot.Average.Niche.squares <- function( annot.table, indexCluster, title )
{
  # initialize the data table 
  grid_data <- data.table(  
    x = c(1,2,3,3,3,2,1,1,2),
    y = c(3,3,3,2,1,1,1,2,2),
    subset = c("red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "brown"),
    proportion = c(0.8, 0.5, 0.7, 0.6, 0.9, 0.4, 0.7, 0.3, 0.8) # Arbitrary proportions for each square
  )
  
  # make the center the index cell 
  res <- showNiche.multiSample(annot.table = annot.table, indexCluster = indexCluster, title=indexCluster, plots = F)
  stableColorsL2.fxn <- stableColorsL2
  
  ave <- res %>% group_by(subset) %>% summarise(mean_proportion = mean(proportion, na.rm = TRUE))
  ave <- ave[1:9,]
  ave$x[1:9] <- grid_data$x; ave$y[1:9] <- grid_data$y
  grid_data <- as.data.table(ave)
  grid_data[x == 2 & y == 2, subset:= indexCluster]
  grid_data[x == 2 & y == 2, mean_proportion:= 0]
  
  # grid_data$mean_proportion <- grid_data$mean_proportion 
  grid_data$label <- gsub(pattern = "_", replacement = " ", x = grid_data$subset)
  grid_data$label <- stringr::str_wrap(grid_data$label, 5)
  grid_data$label.annot.x <- c(0.5, 2, 3.5, 3.5, 3.5, 2, 0.5, 0.5, 2)
  grid_data$label.annot.y <- c(3.75, 3.75, 3.75, 2, 0.25, 0.25, 0.25, 2, 2)
  hex.coord <- data.frame(x = c(2.5, 2.25, 1.75, 1.5, 1.75, 2.25), y=c(2, 2.433, 2.433, 2, 1.567, 1.567))
  print(
    ggplot(grid_data[1:8,], aes( x=x, y=y, fill=subset, label=label, size=mean_proportion)) + 
      geom_polygon(data = hex.coord, inherit.aes = F, aes(x = x, y = y), fill = "white", color = "black", size=0.05) +
      geom_point(shape = 22, color = 'gray50') +   
      geom_text(aes(x = label.annot.x, y=label.annot.y, label = label), size=2) + 
      coord_equal() + # Ensure squares are square (base aspect ratio)
      scale_size(range = c(0,20)) + scale_fill_manual(values = stableColorsL2.fxn) + 
      ggtitle(title) + guides(fill = "none")  +#  scale_size_area() + 
      annotate("text", x = 2, y = 2, label = title, colour = "black", fontface = 'bold', size=2) +  
      theme_void() +
      xlim(c(-0.5,4.5)) + ylim(c(-0.5,4.5))
  )
}

plot.Average.Niche.igraph <- function( annot.table, indexCluster, title )
{
  library(ggraph)
  # make the center the index cell 
  res <- showNiche.multiSample(annot.table = annot.table, indexCluster = indexCluster, title=indexCluster, plots = F)
  stableColorsL2.fxn <- stableColorsL2
  
  ave <- res %>% group_by(subset) %>% summarise(mean_proportion = mean(proportion, na.rm = TRUE))
  ave$subset <- as.character(ave$subset)
  if (indexCluster %in% ave$subset)   {    
    ave$subset[which(ave$subset == indexCluster)] <- paste0(indexCluster, ".")  
    duplColor <- stableColorsL2.fxn[which(names(stableColorsL2.fxn)==indexCluster)]
    names(duplColor) <- paste0(names(duplColor), ".")
    stableColorsL2.fxn <- c(stableColorsL2.fxn, duplColor)
  }
  center_node_df <- data.frame(subset = indexCluster, mean_proportion = 0) 
  setorder(ave, -mean_proportion)  # order ave from biggest to smallest 
  ave <- ave[1:8,]  # keep the top 8 
  grid_data <- as.data.table(bind_rows(center_node_df, ave))  # add the index cell to the data frame 
  grid_data$label <- grid_data$subset  # wrap labels  
  grid_data$label <- gsub(pattern = "_",replacement = " ", x = grid_data$label)
  grid_data$label <- stringr::str_wrap(grid_data$label, 8)  # wrap labels  
  edges <- data.frame(from = indexCluster, to = grid_data$subset[grid_data$subset != indexCluster] )  # new df with directional proportions
  nodes <- data.frame(name = grid_data$subset, label = grid_data$label, size = grid_data$mean_proportion) # map to node size
  graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)  # create an igraph object with nodes and edges 
  print( 
    a <- ggraph(graph, layout = 'dendrogram', circular = TRUE) +
      geom_edge_fan(width = 0.5, color = "gray60", alpha = 0.8) +
      geom_node_point(aes(fill = name, size = size), shape = 21, color = "black", stroke = 0.2) +
      geom_node_text(aes(label = label, x = x * 1.3, y = y * 1.3), repel = F, size = 2, 
                     hjust='outward', vjust='outward') +
      scale_size_continuous( name = "Interaction\nProportion", limits = c(0, 0.5), range = c(0.1, 15), 
                             breaks=c(0, 0.01, 0.05,0.1,0.2,0.4) ) + 
      scale_fill_manual(values = stableColorsL2.fxn, guide = "none") +
      geom_node_label( aes(label = label, filter = name == indexCluster), size = 2, alpha=1, fontface = 'bold') +
      theme_graph(base_family = 'sans') +
      labs(title = title) +
      xlim(c(-3, 3)) + 
      ylim(c(-3, 3)) +
      coord_fixed()
  )
  return(a)
}


findImportantInteractions <- function( interaction.dt, condition.exp, condition.ctrl, ef.threshold, leaveOutTumor = T)
{
  working.dt <- copy(interaction.dt)
  # example output
  #                                                  unified_int type_int original simulations   enrichm p_higher_orig p_lower_orig
  #                                                         <fctr>   <char>    <num>       <num>     <num>         <num>        <num>
  # 1:                                                Tumor--Tumor     homo     5976     110.155  5.748777             0            1
  # 2:    CD8TOXHI_TCF1_HI_TOX_HI_CD8--CD8TOXHI_TCF1_HI_TOX_HI_CD8     homo       61       0.500  5.369234             0            1
  # 3:              TBMIXED_MIXED_GC_TCELL--TBMIXED_MIXED_GC_TCELL     homo     6627     227.565  4.857897             0            1
  #        p.adj_higher p.adj_lower  PI_value int_ranking list_ID condition
  #           <num>       <num>     <num>       <int>  <char>    <char>
  # 1:            0           1  13.22811           1 tumor.1    stage3
  # 2:            0           1  12.35477           2 tumor.1    stage3
  # 3:            0           1  11.17817           3 tumor.1    stage3
  
  # per discussion with Rami Vanguri (6/24), we should consider clipping the enrichm scores at 0. these scores represent deviations away from random but at some level, you cannot be more random than random. the negative enrichm values imply loss of interactions but maybe that has to be thought of distinctly from gain of interactions with respect to the baseline of random interaction 
  working.dt[, enrichm := ifelse(enrichm > 0, yes = enrichm, no=0)]
  delta.pairwise.enrichment<- working.dt %>% 
    group_by(unified_int, condition) %>% 
    summarise(mean = mean(enrichm, na.rm = TRUE), 
              sd = sd(enrichm, na.rm = TRUE), 
              n_replicates = n()) %>%  ungroup()
  # example output
  # A tibble: 1,406 × 5
  #   unified_int                                              condition  mean    sd n_replicates
  #   <fct>                                                    <chr>     <dbl> <dbl>        <int>
  # 1 BCELLS_MANTLE_ZONE_BCELL--TBMIXED_MIXED_GC_TCELL         stage1        0     0            5
  # 2 BCELLS_MANTLE_ZONE_BCELL--TBMIXED_MIXED_GC_TCELL         stage3        0     0           10
  # 3 BCELLS_MIXED_FOLLICULAR_BCELL--Tumor                     stage1        0     0            5
  # 4 BCELLS_MIXED_FOLLICULAR_BCELL--Tumor                     stage3        0     0           10
  # 5 BCELLS_MANTLE_ZONE_BCELL--Tumor                          stage1        0     0            5
  # 6 BCELLS_MANTLE_ZONE_BCELL--Tumor                          stage3        0     0           10
  
  
  # For each group, subtract the 'mean' value of stage1 from stage3. Use `mean[condition == '...']` to select the correct value.
  mean.diff.enrichm <- delta.pairwise.enrichment %>%
    group_by(unified_int) %>%  
    summarise( mean.diff = mean[condition == condition.exp] - mean[condition == (condition.ctrl)], 
               effect.size = (mean[condition == (condition.exp)] - mean[condition == (condition.ctrl)]) / sqrt(0.5*(sd[condition==(condition.exp)]^2 + sd[condition==(condition.ctrl)]^2))) %>% 
    arrange(desc(effect.size))     
  print("mean.diff.enrichm")
  print(as_tibble(mean.diff.enrichm))
  # example output 
  # A tibble: 703 × 3
  #   unified_int                                                            mean.diff       effect.size
  #   <fct>                                                                      <dbl>             <dbl>
  # 1 CD4FOXP3HI_ACTIVATED_TREG--CD8TOXLO_MIXED_TOX_LO_CD8                       0.503              1.61
  # 2 CD4.outlier--CD4FOXP3HI_MIXED_TREG                                         0.885              1.51
  # 3 CD8TOXHI_TCF1_MIXED_TOX_HI_CD8--CD8TOXLO_MIXED_TOX_LO_CD8                  0.499              1.47
  # 4 BCELLS_PC--TBMIXED_BCELL_MYELOID_MIXED                                     0.449              1.47
  # 5 CD8TOXHI_TCF1_MIXED_TOX_HI_CD8--CD8TOXLO_TERMINAL_EFFECTOR_CD8             0.810              1.33
  
  a <- ggplot(mean.diff.enrichm, aes(x = 1:nrow(mean.diff.enrichm), y=effect.size, label=unified_int)) +
    geom_point(size=2, pch=21, color="white",  alpha=0.5, fill="black") +
    ggrepel::geom_text_repel(size=2.5, max.overlaps = 10, segment.color = 'gray60', force=5) +
    ggtitle(paste0("effect size: ", condition.exp, " - ", condition.ctrl))  + theme_classic() +
    theme(axis.text.x = element_blank())
  
  df <- as.data.table(mean.diff.enrichm)[ abs(effect.size) > ef.threshold]
  df <- as.data.table(df)[ !unified_int %like% 'outlier' & !unified_int %like% 'ARTIFACT'] 
  if(leaveOutTumor == T) { df <- as.data.table(df)[!unified_int %like% 'Tumor'] }  # leave out tumor interactions 
  df$unified_int <- factor(df$unified_int, levels = df$unified_int)
  print("effect size table, filtered for effect size threshold")
  print(as_tibble(df))
  # example output
  # # A tibble: 30 × 3
  #    unified_int                                        mean.diff       effect.size
  #    <fct>                                                  <dbl>             <dbl>
  #  1 ACTIVATED_TREG--MIXED_TOX_LO_CD8                       0.503              1.61
  #  2 MIXED_TOX_LO_CD8--TCF1_MIXED_TOX_HI_CD8                0.499              1.47
  #  3 TCF1_MIXED_TOX_HI_CD8--TERMINAL_EFFECTOR_CD8           0.810              1.33
  #  4 CD4_MYELOID_NONIMMUNE_MIXED--TCF1_MIXED_TOX_HI_CD8     0.221              1.31
  #  5 GZMB_CD8--NONIMMUNE_CELL                               0.179              1.18
  #  6 ACTIVATED_MEMORY_CD4--MIXED_TOX_LO_CD8                 0.441              1.18
  #  ...
  # 26 MEMORY_CD8--MIXED_TREG                               -0.505              -1.54
  # 27 ARTIFACT--CD4_MYELOID_NONIMMUNE_MIXED                -0.916              -1.60
  # 28 ARTIFACT--MIXED_TOX_LO_CD8                           -0.530              -1.67
  # 29 CD8_PC_MIXED--GZMB_CD8                               -0.497              -1.68
  # 30 ACTIVATED_TREG--MEMORY_CD8                           -0.778              -2.08
  
  # annotateX.exp <- max(df$effect.size)-0.2*max(df$effect.size)
  # annotateX.ctrl <- min(df$effect.size)+0.2*min(df$effect.size)
  
  # b <- ggplot(df, aes(x = effect.size, y = unified_int))  +   
  #   geom_bar(stat='identity', width=0.05) +  geom_point(size=4, pch=21, color='white', fill='black', stroke=0.2) + 
  #         theme_classic() + ggtitle( paste0("Effect size (filtered): ", condition.exp, " - ", condition.ctrl))  +       
  #   theme(axis.text.y = element_text( size=8), axis.title.y = element_blank(), plot.title=element_text(size=10), 
  #         axis.title.x = element_text(size=10), axis.text.x = element_text(size=10)) + 
  #   theme(legend.position = 'none') + annotate("text", x = annotateX.exp, y = -2.2, size = 5, label = condition.exp) + 
  #   annotate("text", x = annotateX.ctrl, y = -2.2, size = 5, label = condition.ctrl) + xlab("Effect size") + 
  #   coord_cartesian(clip = 'off', ylim= c(0, nrow(df)))
  b <- findImportantInteractions.plotFiltered(df, condition.exp, condition.ctrl)
  return(list(df=df, fullggplot=a, filteredggplot=b)) 
}

findImportantInteractions.plotFiltered <- function(df, condition.exp, condition.ctrl, labelPlot =T)
{
  annotateX.exp <- max(df$effect.size)-0.2*max(df$effect.size)
  annotateX.ctrl <- min(df$effect.size)+0.2*min(df$effect.size)
  
  b <- ggplot(df, aes(x = effect.size, y = unified_int))  +   
    geom_bar(stat='identity', width=0.05, color='grey70', fill='grey50') +  
    geom_point(size=4, pch=21, color='white', fill='grey50', stroke=0.2) + 
    theme_classic() + 
    ggtitle( paste0("Effect size (filtered): ", condition.exp, " vs ", condition.ctrl))  +       
    theme(axis.text.y = element_text( size=8, color='black'), axis.title.y = element_blank(), 
          plot.title=element_text(size=10), 
          axis.title.x = element_text(size=10, color='black'), axis.text.x = element_text(size=10, color='black')) + 
    theme(legend.position = 'none')  
  
  if(labelPlot == T)
    (
      b <- ggplot(df, aes(x = effect.size, y = unified_int))  +   
        geom_bar(stat='identity', width=0.05, color='grey70', fill='grey50') +  
        geom_point(size=4, pch=21, color='white', fill='grey50', stroke=0.2) + 
        theme_classic() + 
        ggtitle( paste0("Effect size (filtered): ", condition.exp, " - ", condition.ctrl))  +       
        theme(axis.text.y = element_text( size=8, color='black'), axis.title.y = element_blank(), 
              plot.title=element_text(size=10), 
              axis.title.x = element_text(size=10, color='black'), axis.text.x = element_text(size=10, color='black')) + 
        theme(legend.position = 'none') + annotate("text", x = annotateX.exp, y = -2.2, size = 5, label = condition.exp) + 
        annotate("text", x = annotateX.ctrl, y = -2.2, size = 5, label = condition.ctrl) + xlab("Effect size")  +
        coord_cartesian(clip = 'off', ylim= c(0, nrow(df)))
    )
  return(b)
}



globalChanges.Effectsize.plot <- function(importIntxn.obj, colNames = names(stableColorsL2), title = "Stage 3 vs Stage 1")
{
  importIntxn.obj <- as.data.table(importIntxn.obj$fullggplot$data)
  res_list <- lapply(X = colNames, 
                     FUN = function(x) { importIntxn.obj[unified_int %like% x, effect.size]} )
  names(res_list) <- colNames
  res_list <- lapply(res_list, as.data.table)
  b <- purrr::list_rbind(res_list, names_to = 'subset')
  names(b)[2] <- "Effect.size"
  b <- b[!is.nan('Effect.size')]; b <- b[!subset %like% 'ARTIFACT']; b<- b[!subset %like% 'outlier']; b <- b[!subset %like% 'Tumor']
  b %>% rstatix::dunn_test(formula = Effect.size ~ subset) 
  order <- b %>% group_by(subset) %>% rstatix::get_summary_stats(type = 'mean') %>% arrange(mean)
  b$subset <- factor(b$subset, levels = order$subset) 
  a <- ggplot(b, aes(y = subset, x=Effect.size)) + geom_vline(xintercept = 0, color='grey50')+ 
    geom_violin() +# geom_point(pch=21, color='grey50',fill='grey90', size=0.5) +  
    theme_bw() + theme(axis.title.y = element_blank()) + theme(axis.text = element_text(color='black')) + 
    ggtitle(title)
  return(a)
}  



plot.boxplot <- function(df, probeSubset, grouping.var, plot.var='proportion', nonparam = T, colorList, title)
{
  df <- df[subset == probeSubset]
  # pVal <- df %>% rstatix::welch_anova_test(formula = get(plot.var) ~ get(grouping.var)) 
  formula.str <- as.formula(paste(plot.var, "~", grouping.var))
  if(nonparam==F)   { pVal <- df %>% rstatix::welch_anova_test(formula = formula.str) }
  if(nonparam==T)   { pVal <- df %>% rstatix::wilcox_test(formula = formula.str)  }
  colorList <- colorList[names(colorList) %like% probeSubset]
  ggplot(df, aes(x = get(grouping.var), y=get(plot.var))) + geom_boxplot(aes(fill = color_key)) + 
    geom_point(size=4, pch=21, fill='white', color='grey50',) + theme_bw() + theme(axis.text = element_text(color='black')) + 
    scale_fill_manual(values = colorList)  + theme(legend.position = 'null', axis.title.x = element_blank()) + 
    annotate('text', x=2, y=max(df$proportion), label = paste("P=", round(pVal$p, 2)), color='black')+
    ylab(paste(probeSubset, plot.var)) + ggtitle(title) 
}




plot.enrichm.interaction <- function(interaction.dt, probeSubset, nonparam=T, plot.var='enrichm', colorList, title)
{
  interaction.dt <- interaction.dt[unified_int == probeSubset]
  formula.str <- as.formula(paste(plot.var, "~", 'condition'))
  if(nonparam==F)   { pVal <- interaction.dt %>% rstatix::welch_anova_test(formula = formula.str) }
  if(nonparam==T)   { pVal <- interaction.dt %>% rstatix::wilcox_test(formula = formula.str)   }
  colorList <- colorList[names(colorList) %like% probeSubset]
  ggplot(interaction.dt, aes(x = condition, y=get(plot.var))) + geom_boxplot(fill = 'grey80', alpha=0.5) + 
    geom_point(size=4, pch=21, fill='white', color='grey30',) + theme_bw() + theme(axis.text = element_text(color='black')) + 
    scale_fill_manual(values = colorList)  + theme(legend.position = 'null', axis.title.x = element_blank()) + 
    annotate('text', x=1, y=1.2*max(interaction.dt[,.(get(plot.var))]), 
             label = paste("P=", round(pVal$p, digits=2)), color='black')+
    ggtitle(title) + ylab("Log2FC enrichment relative to random") #+   scale_y_continuous(expand = expansion(mult = 0.1))
}



neighbors.histogram <- function(target.cell, interacting.type, feature, title, xlab, binwidth=0.25, permutations=10000, 
                                annotate.location = c(12,5000))
{
  target.subset <- target.cell[get(interacting.type) >= 1]  # start with interacting.type neighbors
  val.1 <- Giotto::getExpression(gobject=sln, values='normalized')@exprMat[feature, target.subset$from]    # get specific expression data
  target.subset <- target.cell[get(interacting.type) == 0]  # now with no interacting.type neighbors
  val.0 <- Giotto::getExpression(gobject=sln, values='normalized')@exprMat[feature, target.subset$from]    # get specific expression data
  dt <- data.table(data = c(val.0, val.1), 
                   source = c(rep("No.Neighbors", length(val.0)), 
                              rep("Neighbors", length(val.1))) )  # make a new dt where the expression values of FEAT can be plotted 
  p.val <- permutation.test.neighborhoods(data.tbl = dt, group.0 = "No.Neighbors", group.1 = "Neighbors", num.perm = permutations)
  dt %>% group_by(source) %>% rstatix::get_summary_stats(data, type = 'common')
  num.NoNeighbors <- dt[source == 'No.Neighbors', .N]
  num.Neighbors <- dt[source == 'Neighbors', .N]
  if(num.NoNeighbors > num.Neighbors)        {   dt$source <- factor(dt$source, levels = c("No.Neighbors","Neighbors"))  } 
  else {     dt$source <- factor(dt$source, levels = c("Neighbors","No.Neighbors"))  }
  plot <- ggplot(dt, aes(x = data, fill = source)) + geom_histogram(binwidth = binwidth, alpha=0.6, position='identity') + 
    scale_fill_manual(values = c('grey50','firebrick')) + 
    annotate('text', x=annotate.location[1], y=annotate.location[2], color='black', 
             label= paste0("Mean Neighbors: ", round(mean(dt[source == 'Neighbors', data]), 1), 
                           "\n",
                           "Mean No.Neighbors: ", round(mean(dt[source == 'No.Neighbors', data]), 1),
                           "\nPermutation test:  P <= ", ifelse(p.val==0, 1/permutations, p.val)),
             hjust = 0, size=3) +
    ggtitle(title) + xlim(c(0,28 ))+ 
    xlab(xlab) + ylab("Number of cells") +
    theme_classic() + theme(legend.justification = c("right", "top"), legend.position = 'inside', 
                            legend.title=element_blank()) 
  return(list(plot, dt, p.val))
}







