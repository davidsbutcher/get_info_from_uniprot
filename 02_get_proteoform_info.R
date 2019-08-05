# Excel sheet names are limited to 30 characters.
# This script truncates the filenames after 30
# characters to use as Excel sheet names.

# Packages ------------------------------------------------------------------------------------

library(UniProt.ws)
library(here)
library(glue)
library(svMisc)
library(tidyverse)
library(Peptides)
library(magrittr)
library(writexl)
library(readxl)
library(beepr)
library(RSQLite)
library(DBI)
library(zeallot)
library(GO.db)

# Initial Parameters --------------------------------------------------------------------------

# Specify false discovery rate to use for
# rejection of hits (as decimal) when using a
# tdReport as input - 0.01 is 1% FDR

fdr <- 0.01

# Specify UniProt taxon number to search.
# 83333 -> E. coli K12

UPtaxon <- UniProt.ws(83333)

# Functions -----------------------------------------------------------------------------------

read_tdreport <- function(tdreport, fdr_cutoff = 0.01) {
  
  setwd(here("input"))
  
  #Establish database connection
  
  con <- dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)
  
  # Initialize the tibble where we'll be dropping our most 
  # important data, add NA columns for later use
  
  biologicalproteoform <- con %>%
    RSQLite::dbGetQuery("SELECT Id, ProteoformRecordNum, 
                        IsoformId, ChemicalProteoformId 
                        FROM BiologicalProteoform") %>%
    as_tibble %>% 
    add_column(Qvalue = NA) %>% 
    add_column(filename = NA) %>% 
    add_column(MonoisotopicMass = NA) %>% 
    add_column(AverageMass = NA)
    
  
  #Get accession numbers to correlate with isoform IDs later
  
  accession_nums <- con %>%
    RSQLite::dbGetQuery("SELECT Id, AccessionNumber 
                        FROM Isoform") %>%
    as_tibble
  
  names(accession_nums) <- c("IsoformId", "AccessionNumber")
  
  biologicalproteoform %<>% left_join(accession_nums, by = c("IsoformId")) %>% 
    dplyr::select(AccessionNumber, everything())
  
  # Get Qvals and other info from "GlobalQualitativeConfidence"
  # table
  
  q_vals <- con %>%
    RSQLite::dbGetQuery("SELECT ExternalId, GlobalQvalue, HitId 
                        FROM GlobalQualitativeConfidence") %>%
    as_tibble %>% filter(.$ExternalId != 0)
  
  # Get info from ChemicalProteoform and BiologicalProteoform tables
  # to allow for assignment of masses
  
  bioprot <- con %>%
    RSQLite::dbGetQuery("SELECT Id, ProteoformRecordNum, 
                        IsoformId, ChemicalProteoformId 
                        FROM BiologicalProteoform") %>%
    as_tibble
  
  chemprot <- con %>%
    RSQLite::dbGetQuery("SELECT Id, MonoisotopicMass, 
                        AverageMass 
                        FROM ChemicalProteoform") %>%
    as_tibble
  
  names(chemprot) <- c("ChemicalProteoformId", "MonoisotopicMass", "AverageMass")
  
  mass_tbl <- left_join(bioprot, chemprot)
  
  for (i in seq_along(biologicalproteoform$Id)) {
    
    # Retrieve BiologicalProteoformIDs for which IsoformID
    # is a match to this tibble row value
    
    # bioprotid <- con %>%
    #   RSQLite::dbGetQuery("SELECT Id, IsoformId 
    #                       FROM BiologicalProteoform") %>%
    #   as_tibble %>%
    #   filter(IsoformId == isoform_id$Id[i])
    
    # Get minimum Q value from among all Q values for this biological
    # proteoform ID then assign it to "Qvalue" column
    
    min_q_val <-
      q_vals$GlobalQvalue[q_vals$ExternalId == biologicalproteoform$Id[i]] %>%
      min(na.rm = TRUE)
    
    biologicalproteoform$Qvalue[i] <- min_q_val
    
    # Get HitId from q_vals corresponding to lowest Q value hit
    # to be used to find information on data file name
    
    datafile_hit_id <- 
      q_vals %>%
      filter(GlobalQvalue == min_q_val) %>%
      .$HitId %>% 
      .[1]
    
    # Pull all HitIds and DataFileIds from "Hit" table and
    # filter down to single value containing datafile number
    # for lowest Q-value hit
    
    datafile_id <- con %>%
      RSQLite::dbGetQuery("SELECT Id, DataFileId 
                          FROM Hit") %>%
      as_tibble %>% 
      filter(Id == {{datafile_hit_id}}) %>% 
      .$DataFileId %>% 
      .[1]
    
    # Lookup datafile_id in "Id" column of "Datafile" table
    # to get corresponding file name
    
    biologicalproteoform$filename[i] <- con %>%
      RSQLite::dbGetQuery("SELECT Id, Name 
                          FROM DataFile") %>%
      as_tibble %>% 
      filter(Id == datafile_id) %>% 
      .$Name
   
    # Get masses from mass_tbl and add to appropriate column
    
    biologicalproteoform$MonoisotopicMass[i] <- 
      mass_tbl %>% 
      filter(Id == biologicalproteoform$Id[i]) %>%
      .$MonoisotopicMass
    
    biologicalproteoform$AverageMass[i] <- 
      mass_tbl %>% 
      filter(Id == biologicalproteoform$Id[i]) %>%
      .$AverageMass
     
  }
  
  biologicalproteoform_global <- biologicalproteoform
  
  # Filter isoform_id by FDR cutoff value (default 0.01 or 1%)
  # and place results in "output" tibble
  
  output <- biologicalproteoform %>% 
    filter(Qvalue < fdr_cutoff)
  
  setwd(here())
  
  # Close database connection and return output table
  
  dbDisconnect(con)
  
  return(output)
}

getuniprotinfo <- function(tbl, taxon = NULL, tdreport = TRUE) {
  
  # Find column in the input tibble which has "accession"
  # in it and use it to get info from UniProt
  
  accession_name <- grep("accession", names(tbl),
                         ignore.case = TRUE, value = TRUE)
  
  # If the data is coming from a TDreport, add a filename column
  # and Qvalue column
  
  if (tdreport == TRUE) {
    
    results <- dplyr::select(tbl, c(1, 3, 6, 7, 8, 9))
    names(results) <- c("UNIPROTKB", "PFR", "Qvalue", "filename",
                        "MonoisotopicMass", "AverageMass")
    
  } else {
    results <- tibble(UNIPROTKB = pull(tbl, accession_name))
  }
  
  results_temp <<- tibble()
  
  for (i in seq_along(pull(tbl, accession_name))) {

    # Create a safe version of UniProt.ws which will not crash
    # the whole damn program if it fails
    
    safeselect <- safely(UniProt.ws::select)
    
    # For each accession number, pull relevant info from UniProt 
    # and make a new tibble with just that new line
    
    results_newline <<- 
      safeselect(taxon,
                         keys = pull(tbl, accession_name)[i],
                         columns = c("ENTRY-NAME", "GENES",
                                     "PROTEIN-NAMES", "ORGANISM", 
                                     "ORGANISM-ID", "SEQUENCE", 
                                     "FUNCTION",
                                     "SUBCELLULAR-LOCATIONS", 
                                     "GO-ID"),
                 keytype = "UNIPROTKB") 
    
    # If safeselect result evaluates to NULL, columns corresponding
    # to UniProt data are filled with NA
    
    if (is.null(results_newline[["result"]]) == TRUE) {
      
      
      results_newline[["result"]] <- 
        tibble(UNIPROTKB = pull(tbl, accession_name)[i],
               `ENTRY-NAME` = NA,
               `GENES` = NA,
               `PROTEIN-NAMES` = NA,
               `ORGANISM` = NA, 
               `ORGANISM-ID` = NA,
               `SEQUENCE` = NA, 
               `FUNCTION` = NA,
               `SUBCELLULAR-LOCATIONS` = NA, 
               `GO-ID` = NA)
      
    }
    
    results_temp <-  results_newline[["result"]] %>% 
      as_tibble %>% 
      add_column(PFR = results$PFR[i]) %>% 
      union_all(results_temp, results_newline)  
    
    
  }
  
  results_global <<- results
  results_temp_global <<- results_temp
  
  results_after_join <- left_join(results, results_temp)
  return(results_after_join)
  
}

addmasses <- function(tbl) {
  
  # This function uses the Peptides package to determine
  # average and monoisotopic masses based on protein sequence.
  # These values diverge from those in the TDreport by <0.00005 Da.
  
  mutate(tbl, monoiso_mass = mw(tbl$SEQUENCE, monoisotopic = TRUE),
         ave_mass = mw(tbl$SEQUENCE, monoisotopic = FALSE))
  
}

getGOterms <- function(tbl) {
  
  library(magrittr)
  library(GO.db)
  
  print("Getting GO subcellular locations...")
  
  temptbl <- tibble()
  
  temptbl <- tbl %>% 
    add_column(GO_term = NA) %>%
    add_column(GO_subcell_loc = NA)
  
  for (i in seq_along(tbl$`GO-ID`)) {
    
    temptbl$GO_term[i] <- unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>% 
      trimws() %>% Term() %>% paste(collapse = "; ")
    
    temptbl$GO_subcell_loc[i] <- unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>%
      trimws() %>% Term() %>% .[. %in% go_locs] %>% paste(collapse = "; ")
    
  }
  
  return(temptbl)
}

getlocations <- function(resultslist) {
  
  counts <- tibble(filename = names(resultslist),
                   proteoform_count = NA,
                   cytosol_count = NA,
                   membrane_count = NA)
  
  for (i in seq_along(resultslist)) {
    
    counts$proteoform_count[i] <- resultslist[[i]] %>% .$UNIPROTKB %>% length()
    
    counts$cytosol_count[i] <- sum(str_detect(resultslist[[i]]$GO_subcell_loc, "cytosol"), na.rm = TRUE) +
      sum(str_detect(resultslist[[i]]$GO_subcell_loc, "cytoplasm"), na.rm = TRUE)
    
    counts$membrane_count[i] <- sum(str_detect(resultslist[[i]]$GO_subcell_loc, "membrane"), na.rm = TRUE)
    
  }
  
  return(counts)
}

# Read Data Files -----------------------------------------------------------------------------

# Files are read from subdirectory called "input". Data files must
# be in csv or xlsx format and have a column of UniProt IDs whose 
# name includes the word "accession" somewhere. Case doesn't matter
# but spelling does.

filelist <- here("input") %>%
  list.files %>% as.list

setwd(here("input"))

extension <- filelist %>% map(tools::file_ext)

if (length(unique(extension)) > 1) {
  
  stop("More than one kind of file. Try again.")
  
} else if (extension[[1]] == "csv") {
  
  print("Reading csv files...")
  proteoformlist  <- filelist %>% map(read_csv)

} else if (extension[[1]] == "xlsx") {
  
  print("Reading xlsx files...")
  proteoformlist  <- filelist %>% map(read_xlsx)

} else if (extension[[1]] == "tdReport") {
  
  print("Reading tdReport...")
  proteoformlist <- filelist %>% map(read_tdreport, fdr_cutoff = fdr)
  tdreport_file <- TRUE
  
} else {
  print("What are you doing with your life?")
  stop()
}

names(proteoformlist) <- filelist

setwd(here())

# QuickGO_annotations_20190708.tsv contains all GO IDs and corresponding
# GO terms associated with cellular components, i.e. subcellular localization.
# Loading it provides a lookup table we can use to get subcell loc. for
# any relevant GO IDs

go_locs <- read_tsv("QuickGO_annotations_20190708.tsv") %>%
  .["GO NAME"] %>% unique() %>% pull()

# Access UniProt ------------------------------------------------------------------------------

print("Getting info from UniProt...")

# keytypes <- UniProt.ws::keytypes(UPtaxon) %>% enframe
# columns <- UniProt.ws::columns(UPtaxon) %>% enframe

results_proteoform <- proteoformlist %>% 
  map(getuniprotinfo, taxon = UPtaxon, tdreport = tdreport_file) %>%
  map(getGOterms)

if (tdreport_file == FALSE) {
  
  print("Adding masses according to UniProt sequences...")
  results_proteoform %<>% map(addmasses)
  
} else {
  print("Skipping addition of average and monoisotopic masses!")
}

results_proteoform[[length(results_proteoform)+1]] <- getlocations(results_proteoform)


# Output --------------------------------------------------------------------------------------

# An xlsx file with sheets corresponding to files in /data is written to
# the project directory

systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
resultsname <- glue("{systime}_proteoform_results.xlsx")

names(results_proteoform) <- unlist(filelist)

for (i in seq_along(names(results_proteoform))) {
  
  names(results_proteoform)[i] %<>% stringr::str_trunc(28) %>% paste(i, "_", ., sep = "")
  
}

setwd(here("output"))

results_proteoform %>%
  writexl::write_xlsx(path = resultsname)

setwd(here())
