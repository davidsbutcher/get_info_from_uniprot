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
library(openxlsx)
library(DBI)
library(RSQLite)
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

# Add a directory with files to scan. By default all subdirectories
# will be checked UNLESS the directory has "deprecated" in its name.
# MAKE SURE TO ADD THE FINAL FORWARD SLASH for directories
# All input files must be csv, xlsx, or tdReport files.
# You can also add the full path to a single file (including extension).

filedir <- 
  c("Z:/ICR/David Butcher/TDReports/201906_EcoliMG1655_GELFrEE/20190612_GELFrEE_F01-08_1run.tdReport")

# Specify false discovery rate to use for rejection of hits (as decimal)
# when using a tdReport as input - 0.01 is 1% FDR

fdr <- 0.01

# Specify UniProt taxon number to search.
# 83333 -> E. coli K12
# 9606 -> Homo sapiens

taxon_number <- 83333

# QuickGO_annotations_20190708.tsv contains all GO IDs and corresponding
# GO terms associated with cellular components, i.e. subcellular localization.
# Loading it provides a lookup table we can use to get subcell loc. for
# any relevant GO IDs

go_locs_file <- "QuickGO_annotations_20190708.tsv"

# Need to run this command for furrr. If 10 is too many sessions for your
# system try 5. If running <10 files, change workers to be equal to number
# of files. workers = 1 is equivalent to not using furrr at all

plan(multisession(workers = 1))

# Should PushBullet be used to notify when the script is finished?

use_PB <- FALSE

# Load data ---------------------------------------------------------------

# Check for a file in the input folder which contains the UniProt.ws
# taxon data. If it exists, load it. Otherwise, download it SAFELY
# trying up to 10 times

if (file.exists(glue("input/UPtaxon{taxon_number}.rds"))) {
  
  message("Found corresponding RDS file in /input. Loading...")
  UPtaxon <- readRDS(glue("input/UPtaxon{taxon_number}.rds"))
  
} else {
  
  # Establish connection to UniProt WS, safely!
  
  safe_UniProt.ws <- safely(UniProt.ws)
  
  message("Trying to connect to UniProt web service...")
  
  safeUP <- safe_UniProt.ws(taxId = taxon_number)
  
  if (is.null(safeUP[["result"]]) == TRUE) message("Connection failed, trying again!")
  
  iteration_num <- 1
  
  while (is.null(safeUP[["result"]]) == TRUE & iteration_num < 11) {
    
    iteration_num <- iteration_num + 1
    
    message(glue("\nTrying to establish database connection, attempt {iteration_num}"))
    safeUP <- safe_UniProt.ws(taxId = taxon_number)
    
  }
  
  UPtaxon <- safeUP[["result"]]
  
}

message("CONNECTION SUCCEEDED.")

# Load file containing locations corresponding to
# GO terms

go_locs <- go_locs_file %>%
  read_tsv() %>%
  .["GO NAME"] %>%
  unique() %>%
  pull()

# Run Scripts ---------------------------------------------------------------------------------

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

# Optional line used to contact any Pushbullet enabled device.
# View ?pbSetup for help

if (use_PB == TRUE) {

pbPost("note", "R Analysis Finished",
       glue("Elapsed time: {totaltime} min \nDir: {filedir}"))
  
}

# Session info for every run is saved to a txt file in the 
# output directory, in case 

sessionInfo() %>%
  capture.output %>%
  writeLines(glue("output/{systime}_sessionInfo.txt"))

# Close the extra R sessions used by future/furrr

plan(multisession(workers = 1))
