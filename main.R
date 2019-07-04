setwd("C:/Users/dkeo/surfdrive/pd_imaging_scn/pd_scn")

# AHBA data directory
ahba_dir <-"C:/Users/dkeo/surfdrive/AHBA_Arlin"
probeInfo <- read.csv(paste0(ahba_dir, "/probe_info_2018-11-18.csv"))

# Source directory with functions
fun_dir <- paste0(getwd(), "/functions")
R.utils::sourceDirectory(fun_dir)

# Useful variables
donorNames <- c("donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697")
names(donorNames) <- donorNames

source("samples_in_networks.R")


Sys.time()
source("pd_scn/ttest.R")
Sys.time()
source("pd_scn/Rdavid.R")
Sys.time()