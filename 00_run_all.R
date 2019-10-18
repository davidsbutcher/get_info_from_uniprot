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
library(tibble)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(tidyr)

# Initialize Parameters -----------------------------------------------------------------------

setwd(here())

# Add a directory with files to scan. By default all subdirectories
# will be checked UNLESS the directory has "deprecated" in its name.
# MAKE SURE TO ADD THE FINAL FORWARD SLASH for directories
# All input files must be csv, xlsx, or tdReport files.
# You can also add the full path to a single file (including extension).

filedir <- 
  c("Z:/ICR/David Butcher/TDReports/EcoliMG1655/
  20190916_EcoliMG1655_GELFrEE_red_alk/
    20190916_EcoliMG1655WCL_M9-O+L-20190515_2runs_CAMSEARCH.tdReport",
    "Z:/ICR/David Butcher/TDReports/EcoliMG1655/
  201909_EcoliMG1655_PEPPI_M9/
  20190907_EcoliMG1655_PEPPI_M9_F01-F09_CAMsearch.tdReport") %>% 
  str_replace_all("\\n *", "")

# Specify false discovery rate to use for rejection of hits (as decimal)
# when using a tdReport as input - 0.01 is 1% FDR

fdr <- 0.01

# Specify UniProt taxon number to search.
# 83333 -> E. coli K12
# 9606 -> Homo sapiens

taxon_number <- 83333

# Should PushBullet be used to notify when the script is finished?

use_PB <- TRUE

# Should a summary be generated for this analysis?

make_report <- TRUE

# Functions ---------------------------------------------------------------

kickout <- function(list) {
  
  # This function removes any element from the list of input files
  # (from root/input) which does not have one of the allowed
  # extensions or which has "deprecated"
  
  allowed_ext <- c("tdReport", "csv", "xlsx")
  
  for (i in rev(seq_along(list))) {
    
    if (!(tools::file_ext(list[[i]]) %in% allowed_ext)) {
      
      list[[i]] <- NULL 
      
    } else if (str_detect(list[[i]], fixed("deprecated", TRUE)) == TRUE) {
      
      list[[i]] <- NULL 
      
    }
  }
  
  return(list)
}

# Load data ---------------------------------------------------------------

# Data files must
# be in csv, xlsx, or tdReport format and have a column of UniProt IDs whose 
# name includes the word "accession" somewhere. Case doesn't matter
# but spelling does.

if (file.exists(filedir) == TRUE) {
  
  filelist <- filedir %>% as.list() %>% kickout
  filedir <- filedir %>% dirname()
  
  names(filelist) <- seq(1, length(filelist))
  
} else {
  
  filelist <- filedir %>%
    list.files(recursive = T, include.dirs = T, full.names = T) %>%
    as.list %>% kickout
  
  names(filelist) <- seq(1, length(filelist))
  
}

# Need to run this command for furrr. If 10 is too many sessions for your
# system try 5. If running <10 files, change workers to be equal to number
# of files. workers = 1 is equivalent to not using furrr at all

if (length(filelist) >= 8) {
  
  plan(multisession(workers = 8))
  
} else {
  
  plan(multisession(workers = length(filelist)))
       
}

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

go_locs <- readRDS("input/GO_subcellular_locations.rds")

# Run Scripts ---------------------------------------------------------------------------------

# Save start time to variable for use in output filenames

systime <- format(Sys.time(), "%Y%m%d_%H%M%S")

tic()

source("01_get_protein_info.R")
source("02_get_proteoform_info.R")
source("03_output_plots.R")

if (make_report == TRUE) {
  
  rmarkdown::render("04_generate_report.Rmd",
                    output_file = glue("output/report/{systime}_report.html"))
  
}

totaltime <- 
  capture.output(toc()) %>%
  str_extract("[0-9]+") %>%
  as.numeric %>%
  `/`(60) %>%
  round(digits = 2)

print(glue("Elapsed time: {totaltime} min"))

# Optional line used to contact any Pushbullet enabled device.
# View ?pbSetup for help

if (use_PB == TRUE) {

pbPost("note", "R Analysis Finished",
       glue("Elapsed time: {totaltime} min \nDir: {filedir}"))
  
}

# Session info for every run is saved to a txt file in the 
# output directory, in case 

if(dir.exists("output/session_info") == FALSE) dir.create("output/session_info")

sessioninfo::session_info() %>%
  capture.output %>%
  writeLines(glue("output/session_info/{systime}_sessionInfo.txt"))

# Close the extra R sessions used by future/furrr

plan(multisession(workers = 1))
