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

# Functions -----------------------------------------------------------------------------------

read_tdreport <- function(tdreport, fdr_cutoff = 0.01) {
  
  setwd(here("input"))
  
  #Establish database connection
  
  con <- dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)
  
  # Initialize the tibble where we'll be dropping our most 
  # important data, add NA columns for later use
    
  isoform_id <- con %>%
    RSQLite::dbGetQuery("SELECT Id, AccessionNumber 
                        FROM Isoform") %>%
    as_tibble %>% 
    add_column(Qvalue = NA) %>% 
    add_column(filename = NA)
  
  # Get Qvals and other info from "LocalQualitativeConfidence"
  # table
  
  q_vals <- con %>%
    RSQLite::dbGetQuery("SELECT ExternalId, Qvalue, HitId 
                        FROM LocalQualitativeConfidence") %>%
    as_tibble
  
  for (i in seq_along(isoform_id$Id)) {
    
    # Retrieve BiologicalProteoformIDs for which IsoformID
    # is a match to this tibble row value
    
    bioprotid <- con %>%
      RSQLite::dbGetQuery("SELECT Id, IsoformId 
                          FROM BiologicalProteoform") %>%
      as_tibble %>%
      filter(IsoformId == isoform_id$Id[i])
    
    # Get minimum Q value from among all Q values for this Isoform
    # then assign it to "Qvalue" column
    
    min_q_val <- 
      q_vals$Qvalue[q_vals$ExternalId %in% bioprotid$Id]
      min(na.rm = TRUE)
    
    isoform_id$Qvalue[i] <- min_q_val
    
    # Get HitId from q_vals corresponding to lowest Q value hit
    # to be used to find information on data file name
    
    datafile_hit_id <- 
      q_vals %>%
      filter(Qvalue == min_q_val) %>%
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
    
    isoform_id$filename[i] <- con %>%
      RSQLite::dbGetQuery("SELECT Id, Name 
                          FROM DataFile") %>%
      as_tibble %>% 
      filter(Id == datafile_id) %>% 
      .$Name
    
  }
  
  # Filter isoform_id by FDR cutoff value (default 0.01 or 1%)
  # and place results in "output" tibble
  
  output <- isoform_id %>% 
    filter(Qvalue < fdr_cutoff)
  
  setwd(here())
  
  # Close database connection and return output table
  
  dbDisconnect(con)
  
  return(output)
}

getuniprotinfo <- function(tbl, taxon = NULL) {
  
  accession_name <- grep("accession", names(proteinlist[[1]]),
                         ignore.case = TRUE, value = TRUE)
  
  # results <- tibble()
  
  results <- tibble(UNIPROTKB = pull(tbl, accession_name),
                    filename = pull(tbl, filename))
  results_temp <- tibble()
  
  for (i in seq_along(pull(tbl, accession_name))) {
    
    results_temp <- 
      UniProt.ws::select(taxon,
                         keys = pull(tbl, accession_name)[i],
                         columns = c("ENTRY-NAME", "PROTEIN-NAMES",
                                     "ORGANISM", "ORGANISM-ID",
                                     "SEQUENCE", "FUNCTION",
                                     "SUBCELLULAR-LOCATIONS", 
                                     "GO-ID"),
                         keytype = "UNIPROTKB") %>% 
      as_tibble %>% 
      union_all(results_temp, .)
    
  }
  
  # results <- tibble(UNIPROTKB = pull(proteinlist[[1]], AccessionNumber),
  #                   filename = pull(proteinlist[[1]], filename))
  # 
  # results_temp_global <<- results_temp
  
results_after_join <<- left_join(results, results_temp)
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
  proteinlist  <- filelist %>% map(read_csv)
  names(proteinlist) %<>% stringr::str_trunc(30)
  
} else if (extension[[1]] == "xlsx") {
  
  print("Reading xlsx files...")
  proteinlist  <- filelist %>% map(read_xlsx)
  names(proteinlist) %<>% stringr::str_trunc(30)
  
} else if (extension[[1]] == "tdReport") {
  
  print("Reading tdReport...")
  
  # The read_tdreport function is deprecated for now. Instead I have
  # to run this huge code block in the global environment!!!
  # 
  # proteinlist <- filelist %>% map(read_tdreport, fdr_cutoff = fdr)
  
  proteinlist <- list()
  
  for (i in seq_along(filelist)) {
    
    setwd(here("input"))
    
    #Establish database connection
    
    con <- dbConnect(RSQLite::SQLite(), ":memory:", dbname = filelist[[i]])
    
    # Initialize the tibble where we'll be dropping our most 
    # important data, add NA columns for later use
    
    isoform_id <- con %>%
      RSQLite::dbGetQuery("SELECT Id, AccessionNumber 
                        FROM Isoform") %>%
      as_tibble %>% 
      add_column(Qvalue = NA) %>% 
      add_column(filename = NA)
    
    # Get Qvals and other info from "LocalQualitativeConfidence"
    # table
    
    q_vals <- con %>%
      RSQLite::dbGetQuery("SELECT ExternalId, Qvalue, HitId 
                        FROM LocalQualitativeConfidence") %>%
      as_tibble
    
    for (j in seq_along(isoform_id$Id)) {
      
      # Retrieve BiologicalProteoformIDs for which IsoformID
      # is a match to this tibble row value
      
      bioprotid <- con %>%
        RSQLite::dbGetQuery("SELECT Id, IsoformId 
                          FROM BiologicalProteoform") %>%
        as_tibble %>%
        filter(IsoformId == isoform_id$Id[j])
      
      # Get minimum Q value from among all Q values for this Isoform
      # then assign it to "Qvalue" column
      
      min_q_val <- 
        q_vals$Qvalue[q_vals$ExternalId %in% bioprotid$Id] %>% 
        min(na.rm = TRUE)
      
      isoform_id$Qvalue[j] <- min_q_val
      
      # Get HitId from q_vals corresponding to lowest Q value hit
      # to be used to find information on data file name
      
      datafile_hit_id <- 
        q_vals %>%
        filter(Qvalue == min_q_val) %>%
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
      
      isoform_id$filename[j] <- con %>%
        RSQLite::dbGetQuery("SELECT Id, Name 
                          FROM DataFile") %>%
        as_tibble %>% 
        filter(Id == datafile_id) %>% 
        .$Name
      
    }
    
    # Filter isoform_id by FDR cutoff value (default 0.01 or 1%)
    # and place results in "output" tibble
    
    proteinlist[[i]] <- isoform_id %>% 
      filter(Qvalue < fdr)
    
    setwd(here())
    
    # Close database connection and return output table
    
    dbDisconnect(con)
    
  }
  
  names(proteinlist) %<>% stringr::str_trunc(30)
  
  
} else {
  print("What are you doing with your life?")
  stop()
}

names(proteinlist) <- filelist

setwd(here())

# QuickGO_annotations_20190708.tsv contains all GO IDs and corresponding
# GO terms associated with cellular components, i.e. subcellular localization.
# Loading it provides a lookup table we can use to get subcell loc. for
# any relevant GO IDs

go_locs <- read_tsv("QuickGO_annotations_20190708.tsv") %>%
  .["GO NAME"] %>% unique() %>% pull()

# Access UniProt ------------------------------------------------------------------------------

# Specify UniProt taxon number to search.
# 83333 -> E. coli K12

UPtaxon <- UniProt.ws(83333)

keytypes <- UniProt.ws::keytypes(UPtaxon) %>% enframe
columns <- UniProt.ws::columns(UPtaxon) %>% enframe

results <- proteinlist %>% 
  map(getuniprotinfo, taxon = UPtaxon) %>%
  map(as_tibble) %>%
  map(getGOterms) %>% 
  map(addmasses)


#   count_cytosol <- sum(str_detect(results[[i]]$GO_subcell_loc, "cytosol"), na.rm = TRUE) + 
#     sum(str_detect(results[[i]]$GO_subcell_loc, "cytoplasm"), na.rm = TRUE)
#   count_membrane <- sum(str_detect(results[[i]]$GO_subcell_loc, "membrane"), na.rm = TRUE)


# Output --------------------------------------------------------------------------------------

# An xlsx file with sheets corresponding to files in /data is written to
# the project directory

systime <- format(Sys.time(), "%Y%m%d_%H%M%S")
resultsname <- glue("{systime}_results.xlsx")

names(results) %<>% stringr::str_trunc(30)

results %>%
  writexl::write_xlsx(path = resultsname)

beepr::beep(sound = 8)

# results_plot <- ggplot() +
#   # geom_histogram(data = results[[1]], aes(ave_mass),
#   #                binwidth = 1000, col = "black", fill = "blue", alpha = 0.4) +
#   # geom_histogram(data = results[[2]], aes(ave_mass),
#   #                binwidth = 1000, col = "black", fill = "red", alpha = 0.4) +
#   # geom_histogram(data = results[[3]], aes(ave_mass),
#   #                binwidth = 1000, col = "black", fill = "yellow", alpha = 0.4) +
#   geom_freqpoly(data = results[[1]], aes(ave_mass), color = "red", binwidth = 1000) +
#   geom_freqpoly(data = results[[2]], aes(ave_mass), color = "blue", binwidth = 1000) +
#   geom_freqpoly(data = results[[3]], aes(ave_mass), color = "dark green", binwidth = 1000) +
#   geom_freqpoly(data = results[[4]], aes(ave_mass), color = "black", binwidth = 1000) +
#   theme_minimal()

