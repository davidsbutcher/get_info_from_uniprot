# This script will pull all proteins/proteoforms with Qvalues
# below the provided FDR cutoff (defaults to 0.01) and print the
# results to Excel files.
#
# Input can be tdReport, csv, or xlsx files (single sheet).
# csv and xlsx files must have a column which is a list of
# UniProt accession numbers with a column name that includes
# the word "Accession" SOMEWHERE. Capitalization doesn't
# matter but spelling does. 

library(here)
library(tictoc)
library(glue)
library(rstudioapi)
library(RPushbullet)
library(magrittr)

setwd(here())

tic()

source("01_get_protein_info.R")
source("02_get_proteoform_info.R")
source("03_output_plots.R")

totaltime <- capture.output(toc()) %>%
  str_extract("[0-9]+") %>%
  as.numeric %>% 
  `/`(60) %>% round(digits = 2)

print(glue("Elapsed time: {totaltime} min"))

pbPost("note", "R Analysis Finished",
       glue("Elapsed time: {totaltime} min \nDir: {filedir}")) 
