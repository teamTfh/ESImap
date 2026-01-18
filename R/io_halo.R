

#' Clean HALO Data
#'
#' Parses a raw HALO data object, extracting expression data, metadata, and spatial locations.
#'
#' @param halo.obj A data.frame or data.table read from a HALO csv file.
#' @param label A character string to prefix cell IDs. Defaults to a random 5-letter string.
#' @param saveColumns A character vector of column names to retain as metadata.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{expr}: Matrix of expression values.
#'   \item \code{metadata}: data.table of cell metadata.
#'   \item \code{locations}: data.table of X/Y coordinates.
#' }
#' @import data.table
#' @export
cleanHALO.broken <- function(halo.obj, label=paste0(sample(letters, 5, T), collapse = ""), saveColumns)
{
  # Ensure input is a data.table
  if (!data.table::is.data.table(halo.obj)) {   data.table::setDT(halo.obj)   }

  # Clean column names
  data.table::setnames(halo.obj, names(halo.obj), gsub(" ", ".", names(halo.obj)))

  # Extract Expression Data
  cellIntensity <- grep("Cell.Intensity", names(halo.obj), value = TRUE)
  expr <- setDF(halo.obj[, .SD, .SDcols = cellIntensity])
  names(expr) <- gsub("\\.Cell.Intensity$", "", names(expr))
  rownames(expr) <- paste0(label, "_", halo.obj$Object.Id)

  # Extract Locations
  locations <- data.table::data.table(
    x = rowMeans(halo.obj[, c('XMax', 'XMin'), with = FALSE]),
    y = rowMeans(halo.obj[, c('YMax', 'YMin'), with = FALSE])
  )

  # Extract Metadata
  saveColumns <-   gsub(" ", ".", saveColumns)
  metadata <- halo.obj[, .SD, .SDcols = saveColumns]
  metadata[, cellID := paste0(label, "_", Object.Id)]
  metadata[, x := locations$x]
  metadata[, y := locations$y]

  # Clean NaNs
  for (j in seq_len(ncol(metadata))) {
    set(metadata, which(is.nan(metadata[[j]])), j, 0)
  }

  return(list(expr = t(expr), metadata = metadata, locations = locations))
}




#' Create and Normalize Giotto Object from HALO Data
#'
#' Wraps the cleaning, creation, filtering, and normalization steps into one pipeline.
#'
#' @param inputFile A raw data.table or data.frame from HALO.
#' @param label A character string for the sample name.
#' @param instructions A list of Giotto instructions (created via \code{Giotto::createGiottoInstructions}).
#'
#' @return A processed Giotto object.
#' @import Giotto
#' @export
makeGiottoObject.broken <- function(inputData, label, instructions=NULL)
{
  # 1. Clean the HALO data
  out <- inputData
  message(paste0("Making the Giotto object for ", label))

  # 2. Create Object
  # Note: using Giotto:: prefix to ensure the package finds the function
  sln <- Giotto::createGiottoObject(expression = out$expr,spatial_locs = out$locations,instructions = instructions)

  message("Adding cell metadata")
  sln <- Giotto::addCellMetadata(sln, new_metadata = out$metadata, by_column = TRUE, column_cell_ID = "cellID")

  # 3. Filter
  cell_metadata <- Giotto::pDataDT(sln)
  keepID <- cell_metadata[Classifier.Label != "Background" & Classifier.Label != "Artifact", ]$cell_ID
  sln <- Giotto::subsetGiotto(sln, cell_ids = keepID)

  # 4. Normalize
  message("Starting normalization")
  sln <- Giotto::normalizeGiotto(gobject = sln, scalefactor = 6000, norm_methods = "standard", verbose = TRUE, log_norm = FALSE,
    library_size_norm = FALSE, scale_feats = FALSE, scale_cells = TRUE)
  message("Finished normalization")

  sln <- Giotto::addStatistics(gobject = sln, expression_values = "normalized")

  return(sln)
}


