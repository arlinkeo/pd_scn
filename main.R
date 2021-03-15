# Main script to run all analyses

setwd("C:/Users/Arlin/surfdrive/pd_imaging_scn/pd_scn")
options(stringsAsFactors = FALSE)

# Useful variables
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames
network_names <- c(Network_C = "Posterior cingulate network", 
                   Network_D = "Anterior cingulate network")

# Source directory with functions
fun_dir <- paste0(getwd(), "/functions/")
R.utils::sourceDirectory(fun_dir, modifiedOnly = FALSE)

# Make output folder
dir.create("output")

# AHBA Data preprocessing
source("pd_braak/probe2gene_AHBA.R")

# Select AHBA samples within SCNs
source("samples_in_networks.R")

# Differential expression between networks
source("differential_expression.R")

# Differential stability of DEGs
source("differential_stability.R")

# Gene set enrichment analysis
source("pathway_analysis.R")
source("go_enrichment.R")

# Cell-type marker enrichment
source("celltype_marker_conversion.R")
source("celltype_enrichment_degs.R")
source("heatmap_celltypes.R")
source("ewce.R")

# Enrichment of disease-associated genes
source("disease_enrichment_degs.R")