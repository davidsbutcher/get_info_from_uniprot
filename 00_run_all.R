library(here)
library(tictoc)
library(glue)
library(rstudioapi)
library(RPushbullet)

setwd(here())

tic()
source("01_get_protein_info.R")
source("02_get_proteoform_info.R")
source("03_output_plots.R")
tictoc_time <- capture.output(toc())

pbPost("note", "Analysis Finished", tictoc_time) 
