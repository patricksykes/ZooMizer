## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and three fish groups
## This code is to run the model across multiple cores
##
## Updated for ZooMizer - 18th August 2020

library(mizer) 

source("fZooMizer_run.R") #source the model code

environment(new_project_simple) <- asNamespace('mizer')
assignInNamespace("project_simple", new_project_simple, ns = "mizer")

jobname <- '20200818_ZooMizerTest' #job name used on queue

## Build environmental data
enviro_data <- readRDS("data/envirofull_20200317.RDS")
enviro_data$tmaxx <- 1000

# NOW WE BREAK THIS UP INTO CHUNKS OF 10 RUNS AND CREATE A LIST
n_runs <- 10
enviro_list <- split(enviro_data, ((as.numeric(rownames(enviro_data)) - 1) %/% n_runs) + 1)

HPC <- FALSE # Is this being run on a HPC or will we choose the row

enviro_row <- 1

Groups <- read.csv("data/TestGroups_mizer.csv", stringsAsFactors = FALSE) # Load in functional group information


if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
} else {
  ID <- enviro_row
}

all_params <- enviro_list[[ID]]

for(r in 1:n_runs){
  input_params <- all_params[r,]
  ID_char <- sprintf("%06d",input_params$cellID) # Set the ID as a 4 digit character so it will sort properly
  
  out@params@other_params$runtime <- system.time(
    out <- fZooMizer_run(groups = Groups, input = input_params)
  )
  
  saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))
  rm(out)
}
