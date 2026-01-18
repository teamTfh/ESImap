# ESImap

**Spatial Analysis and Neighborhood Mapping for HALO Data**

ESImap is an R package designed to provide a specialized toolset for processing HALO imaging data using [Giotto](https://github.com/drieslab/Giotto). It enables advanced workflows for data integration, niche composition analysis, and neighborhood-dependent feature expression testing.

## Overview

High-plex imaging platforms like HALO generate complex spatial data. ESImap bridges the gap between raw HALO outputs and advanced spatial transcriptomics/proteomics analysis.

**Key capabilities include:**
*   **Streamlined Data Import:** Parse raw HALO CSVs and convert them directly into Giotto objects with metadata and spatial locations preserved.
*   **Neighborhood Enrichment Analysis:** Compute cell-cell interaction enrichment using permutation-based simulations (`cellProx` workflow).
*   **Niche Visualization:** Visualize cellular neighborhoods using lollipop charts, stacked bars, and network graphs.
*   **Differential Expression in Niches:** Test if feature expression (genes/proteins) changes based on a cell's neighbors.

## Installation

You can install the development version of ESImap from GitHub:

```r
# install.packages("devtools")
devtools::install_github("teamTfh/ESImap")
```

## Dependencies

*   [Giotto](https://rubd.github.io/Giotto_site/)
*   data.table
*   ggplot2
*   dplyr, tidyr, stringr
*   ggraph, tidygraph, ggrepel
*   gridExtra

## Usage

### 1. Data Import & Preprocessing

ESImap provides `cleanHALO` to parse raw data and `makeGiottoObject` to seamlessly create a Giotto object ready for analysis.

```r
library(ESImap)
library(Giotto)

# Load your raw HALO csv file
# halo_df <- fread("path/to/halo_data.csv")

# Create a Giotto object (performs cleaning, filtering, and normalization)
# You can provide Giotto instructions if needed
my_giotto_obj <- makeGiottoObject(inputFile = halo_df,
                                  label = "Sample_01")
```

### 2. Neighborhood Enrichment Analysis

Analyze which cell types interact more frequently than expected by chance using `cellProx.sim` and `cellProx.calcP`.

```r
# Run simulations to generate background distribution of interactions
# Requires a spatial network (e.g., 'spatialknn.k5') to be already computed in the Giotto object
sim_results <- cellProx.sim(gobject = my_giotto_obj,
                            cluster_column = "cell_type",
                            spatial_network_name = "spatialknn.k5",
                            number_of_simulations = 1000)

# Calculate P-values and Enrichment scores
results <- cellProx.calcP(sim_results, number_of_simulations = 1000)

# View enrichment results
head(results$enrichm_res)
```

### 3. Niche Visualization

Visualize the composition of cell neighborhoods.

#### Single Sample Niche
Show the neighbors of a specific cell type (e.g., "Tumor") in one sample.

```r
# Assume 'annot_table' is your annotated spatial network
# showNicheSingleSample(annot.table = annot_table,
#                       indexCluster = "Tumor",
#                       plot_type = "lollipop")
```

#### Network Graph
Visualizes the top neighbors of a specific cell type as a network graph.

```r
# plotAverageNicheIgraph(annot.list = list_of_annotated_networks,
#                        indexCluster = "Macrophage",
#                        top_n = 8)
```

### 4. Neighborhood-Dependent Expression

Test if a marker (e.g., "PDL1") is differentially expressed in cells that have "T-cells" as neighbors versus those that do not.

```r
# Calculate neighborhood proportions first
# neighborhoods <- defineNeighborhoods(annot.list = list_of_annotated_networks)

# Plot histogram and calculate p-value
# plotNeighborHistogram(gobject = my_giotto_obj,
#                       neighborhood_prop_table = neighborhoods$prop.table[[1]],
#                       target_cluster = "Tumor",
#                       neighbor_cluster = "CD8_T_cell",
#                       feature = "PDL1")
```

## Authors

*   **Ramin Herati** - *ramin.herati@nyulangone.org*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
