## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and three fish groups
## This code is to run the model across multiple cores
##
## Updated for this special case - 26th May 2020

# library(Rcpp) # Only needed if we are running with Rcpp code.

source("fZooMSS_Model.R") #source the model code
source("fZooMSS_CalculatePhytoParam.R")

jobname <- '202000930_zoomizercomparison' #job name used on queue

## Build environmental data
enviro <- readRDS("data/enviro_test20.RDS")
enviro_data <- fZooMSS_CalculatePhytoParam(enviro) # Calculate Phytoplankton Parameters
enviro_data$dt
saveRDS(enviro_data, "enviro_Matrix.RDS")

# NOW WE BREAK THIS UP INTO CHUNKS OF 10 RUNS AND CREATE A LIST
n_runs <- 1
enviro_list <- split(enviro_data, ((as.numeric(rownames(enviro_data)) - 1) %/% n_runs) + 1)

HPC <- FALSE # Is this being run on a HPC or will we choose the row
SaveTimeSteps <- TRUE # Should we save all time steps

Groups <- read.csv("TestGroups.csv", stringsAsFactors = FALSE) # Load in functional group information

### No need to change anything below here.
if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
} else {
  ID <- enviro_row
}

all_params <- enviro_list[[ID]]
r(r in 1:n_runs){
  input_params <- enviro_list[[1]]
  
  ID_char <- sprintf("%06d",input_params$cellID) # Set the ID as a 4 digit character so it will sort properly

  out$model$model_runtime <- system.time(
    out <- fZooMSS_Model(input_params, Groups, SaveTimeSteps)
  )

  saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))
  rm(out)
}
