


cleanHALO <- function(halo.obj, label=paste0(sample(letters, 5, T), collapse = ""))
{  # load in a halo object, make a unique column via a random string of 5 letters if no label supplied
  # halo.obj <- inputFile
  expr <- data.frame();   metadata <- data.table(); locations <- data.table()
  setnames(halo.obj, names(halo.obj), gsub(" ", ".", names(halo.obj)))
  cellIntensity <- grep("Cell.Intensity", names(halo.obj), value=T)
  expr <- halo.obj[, which(colnames(halo.obj) %in% cellIntensity)]  # save the Cell.Intensity columns as expression
  names(expr) <- gsub("\\.Cell.Intensity$","",names(expr))  # remove the words ".Cell.Intensity" from columns
  rownames(expr) <- paste0(label, "_", halo.obj$Object.Id)
  locations$x <- rowMeans(halo.obj[,c('XMax', 'XMin')]);   locations$y <- rowMeans(halo.obj[,c('YMax', 'YMin')])
  
  saveColumns <- c(colnames(halo.obj)[grep('Object.Id', colnames(halo.obj), value=F)],
                   colnames(halo.obj)[grep('Classifier.Label', colnames(halo.obj), value=F)],
                   colnames(halo.obj)[grep("Positive.Classification$", colnames(halo.obj), value=F)], 
                   colnames(halo.obj)[9:17])
  metadata <- halo.obj[, saveColumns]
  metadata$cellID <- paste0(label, "_", halo.obj$Object.Id)
  metadata$cellArea <- halo.obj$Cell.Area..µm..
  metadata$x <- locations$x;   metadata$y <- locations$y
  
  metadata[metadata == "NaN"] <- 0
  return(out <- list(expr = t(expr), metadata = metadata, locations = locations))
}


makeGiottoObject <- function(inputFile, label=sample.name)
{
  out <- cleanHALO(inputFile, label)
  print(paste0("Making the Giotto object for ", label))
  sln <- createGiottoObject(expression = out$expr,  spatial_locs = out$locations, instructions = instructions)
  print("Adding cell metadata"); # gc()  
  sln <- addCellMetadata(sln, new_metadata = out$metadata, by_column = TRUE, column_cell_ID = "cellID" )
  cell_metadata <- pDataDT(sln) # pull metadata out into new object
  
  #filter out undesirable cells
  keepID <- cell_metadata[Classifier.Label != "Background" & Classifier.Label != "Artifact", ]$cell_ID
  sln <- subsetGiotto(sln, cell_ids = keepID)
  
  print("Starting normalization")
  sln <- normalizeGiotto(gobject = sln, scalefactor = 6000, norm_methods = "standard",
                         verbose = TRUE, log_norm = FALSE, library_size_norm = FALSE,
                         scale_genes = FALSE, scale_cells = TRUE)
  print("Finished normalization")
  sln <- addStatistics(gobject = sln, expression_values = "normalized")
  return(sln)
}

