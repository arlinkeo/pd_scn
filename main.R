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




# AHBA data directory and data
ahba_dir <-"C:/Users/Arlin/surfdrive/AHBA_Arlin"
probeInfo <- read.csv(paste0(ahba_dir, "/probe_info_2018-11-18.csv"))
brainExpr <- readRDS(paste0(ahba_dir, "/gene_expr.RDS"))
ontology <- read.csv(paste0(ahba_dir, "/Ontology.csv"))
# sample_annot <- lapply(donorNames, function(d){ # Sample info per donor
#   read.csv(paste0("../AHBA_Arlin/sample_info_", d, "_2018-11-18.csv"))
# })



# Run scripts
source("samples_in_networks.R")
source("differential_expression.R")
source("pathway_analysis.R")
source("celltype_marker_conversion.R")
source("celltype_enrichment_degs.R")
source("heatmap_celltypes.R")
source("disease_enrichment_degs.R")