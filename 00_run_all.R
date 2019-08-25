# This script will pull all proteins/proteoforms with Qvalues
# below the provided FDR cutoff (defaults to 0.01) and print the
# results to Excel files.
#
# Input can be tdReport, csv, or xlsx files (single sheet).
# csv and xlsx files must have a column which is a list of
# UniProt accession numbers with a column name that includes
# the word "Accession" SOMEWHERE. Capitalization doesn't
# matter but spelling does. 

# Packages ------------------------------------------------------------------------------------

library(here)
library(svMisc)
library(furrr)
library(Peptides)
library(magrittr)
library(writexl)
library(readxl)
library(RSQLite)
library(DBI)
library(GO.db)
library(tools)
library(progress)
library(tictoc)
library(glue)
library(RPushbullet)
library(ggpubr)
library(UniProt.ws)
library(magrittr)
library(tidyverse)

# Initialize Parameters -----------------------------------------------------------------------

setwd(here())

## Add a directory with files to scan. By default all subdirectories 
## will be checked UNLESS the directory has "deprecated" in its name.
## MAKE SURE TO ADD THE FINAL FORWARD SLASH for directories
## All input files must be csv, xlsx, or tdReport files.
## You can also add the full path to a single file (including extension).

filedir <- 
  c("G:/My Drive/R projects/20190630_get_info_from_uniprot/input/peppi_f03_1dtt_2iaa.csv")

# Specify false discovery rate to use for
# rejection of hits (as decimal) when using a
# tdReport as input - 0.01 is 1% FDR

fdr <- 0.01

# Specify UniProt taxon number to search.
# 83333 -> E. coli K12
# 9606 -> Homo sapiens

UPtaxon <- UniProt.ws(83333)

# QuickGO_annotations_20190708.tsv contains all GO IDs and corresponding
# GO terms associated with cellular components, i.e. subcellular localization.
# Loading it provides a lookup table we can use to get subcell loc. for
# any relevant GO IDs

go_locs_file <- "QuickGO_annotations_20190708.tsv"

# Need to run this command for furrr. If 10
# is too many sessions for your system try 5

plan(multisession(workers = 10))


# Run Scripts ---------------------------------------------------------------------------------

# Load file containing locations corresponding to
# GO terms

go_locs <- go_locs_file %>%
  read_tsv() %>%
  .["GO NAME"] %>%
  unique() %>%
  pull()

# Save start time to variable for use in output filenames

systime <- format(Sys.time(), "%Y%m%d_%H%M%S")

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

sessionInfo() %>%
  capture.output %>%
  writeLines(glue("output/{systime}_sessionInfo.txt"))
