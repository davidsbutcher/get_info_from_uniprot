# Excel sheet names are limited to 30 characters.
# This script truncates the filenames after 30
# characters to use as Excel sheet names.

# Run 01_get_protein_info first!!

# Functions -----------------------------------------------------------------------------------

kickout <- function(list) {
  
  # This function removes any element from the list of input files
  # which does not have one of the allowed
  # extensions
  
  allowed_ext <- c("tdReport", "csv", "xlsx")
  
  for (i in rev(seq_along(list))) {
    
    if (!(tools::file_ext(list[[i]]) %in% allowed_ext)) {
      list[[i]] <- NULL 
    }
  }
  
  return(list)
}

read_tdreport <- function(tdreport, proteinlist, fdr_cutoff = 0.01, file_dir = NULL) {
  
  setwd(file_dir)
  
  #Establish database connection. Keep trying until it works!
  
  safe_dbConnect <- safely(dbConnect)
  
  safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:",
                            dbname = tdreport,
                            synchronous = NULL)
  
  iteration_num <- 1
  
  while (is.null(safecon[["result"]]) == TRUE) {
    
    print(glue("\nTrying to establish database connection, attempt {iteration_num}"))
    safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:",
                              dbname = tdreport,
                              synchronous = NULL)
    
  }
  
  con <- safecon[["result"]]
  
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
    
    qValueExists <- any(!is.na(q_vals$GlobalQvalue[q_vals$ExternalId == biologicalproteoform$Id[i]]))
    
    # If there is a Q value, get minimum Q value from among all Q values
    # for this Isoform ID then assign it to "Qvalue" column. Otherwise
    # assign value 1000000 so it is rejected when performing FDR cutoff
    
    if (qValueExists == TRUE) {
      
      min_q_val <-
        q_vals$GlobalQvalue[q_vals$ExternalId == biologicalproteoform$Id[i]] %>%
        min(na.rm = TRUE)
      
      
    } else {
      
      message(glue("NO VALID Q VALUE FOUND, inserting garbage value in row {i}"))
      
      min_q_val <- 1000000
      
    }
    
    biologicalproteoform$Qvalue[i] <- min_q_val
    
    if (qValueExists == TRUE) {
      
      datafile_hit_id <- 
        q_vals %>%
        filter(GlobalQvalue == min_q_val) %>%
        .$HitId %>% 
        .[1]
      
    } else {
      
      datafile_hit_id <- NA
      
    }
    
    # Pull all HitIds and DataFileIds from "Hit" table and
    # filter down to single value containing datafile number
    # for lowest Q-value hit
    
    if (qValueExists == TRUE & is.na(datafile_hit_id) == FALSE) {
      
      datafile_id <- con %>%
        RSQLite::dbGetQuery("SELECT Id, DataFileId 
                          FROM Hit") %>%
        as_tibble %>% 
        filter(Id == {{datafile_hit_id}}) %>% 
        .$DataFileId %>% 
        .[1]
      
    } else {
      
      datafile_id <- NA
      
    }
    
    # Lookup datafile_id in "Id" column of "Datafile" table
    # to get corresponding file name
    
    if (qValueExists == TRUE & is.na(datafile_id) == FALSE) {
      
      biologicalproteoform$filename[i] <- con %>%
        RSQLite::dbGetQuery("SELECT Id, Name 
                          FROM DataFile") %>%
        as_tibble %>% 
        filter(Id == datafile_id) %>% 
        .$Name
      
    } else {
      
      biologicalproteoform$filename[i] <- NA
      
    }
    
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
  
  # Filter isoform_id by FDR cutoff value (default 0.01 or 1%)
  # and place results in "output" tibble
  
  biologicalproteoform %<>% 
    filter(Qvalue < fdr_cutoff)
  
  revseq <- rev(seq_along(biologicalproteoform$AccessionNumber))
  
  for (i in seq_along(biologicalproteoform$AccessionNumber)) {
    
    # message(glue("\n\n\nBegin iteration {i} of loop.."))
    
    accessiontocheck <- biologicalproteoform$AccessionNumber[revseq[i]]
    
    if (accessiontocheck %in% proteinlist$AccessionNumber == TRUE) {
      
      qvaltocheck <- proteinlist %>% filter(AccessionNumber == accessiontocheck) %>% .$Qvalue
      
    } else {
      
      qvaltocheck <- NA 
      
    }
    
    if (is.na(qvaltocheck) == FALSE & qvaltocheck > fdr_cutoff) {
      
      biologicalproteoform %<>% .[-revseq[i],]
      message(glue("Q value above FDR, erasing proteoform {i}..."))
  
    } else if (is.na(qvaltocheck) == TRUE) {
    
      biologicalproteoform %<>% .[-revseq[i],]
      message(glue("Q value is NA, erasing proteoform {i}..."))
      
    } else {
     
      # print(glue("**NOT** erasing proteoform {i}!"))
       
    }
    
    # print (biologicalproteoform)
    # print(glue("End of iteration {i} of loop."))

  }
  
  setwd(here())
  
  # Close database connection and return output table
  
  dbDisconnect(con)
  
  print("Finished checking corresponding protein entries for Q values above cutoff!")
  print(biologicalproteoform)
  
  return(biologicalproteoform)
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
  
  results_temp <- tibble()
  
  # Prepare progress bar
  
  pb <- progress_bar$new(
    format = "  Accessing UniProt [:bar] :percent eta: :eta :spin",
    total = length(pull(tbl, accession_name)),
    clear = FALSE, width= 60)
  
  for (i in seq_along(pull(tbl, accession_name))) {
    
    # Create a safe version of UniProt.ws which will not crash
    # the whole damn program if it fails
    
    safeselect <- safely(UniProt.ws::select)
    
    # For each accession number, pull relevant info from UniProt 
    # and make a new tibble with just that new line
    
    suppressMessages(results_newline <- 
                       safeselect(taxon,
                                  keys = pull(tbl, accession_name)[i],
                                  columns = c("ENTRY-NAME", "GENES",
                                              "PROTEIN-NAMES", "ORGANISM", 
                                              "ORGANISM-ID", "SEQUENCE", 
                                              "FUNCTION",
                                              "SUBCELLULAR-LOCATIONS", 
                                              "GO-ID"),
                                  keytype = "UNIPROTKB"))
    
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
    
    pb$tick()
  }
  
  results_after_join <- left_join(results, results_temp, by = "PFR")
  
  return(results_after_join)
  
}

getuniprotinfo2 <- function(tbl, filelistnum, tdreport = TRUE) {
  
  # Instead of querying UniProt for information we have already gotten
  # in 01_get_protein_info.R, why not bring it in from results_protein?
  #
  # Find column in the input tibble which has "accession"
  # in it and use it to get info from UniProt
  
  accession_name <- grep("accession", names(tbl),
                         ignore.case = TRUE, value = TRUE)
  
  # If the data is coming from a TDreport, select relevant columns 
  # from input to function
  
  if (tdreport == TRUE) {
    
    results <- dplyr::select(tbl, c(1, 3, 6, 7, 8, 9))
    names(results) <- c("UNIPROTKB", "ProteoformRecordNum", "Qvalue", "filename",
                        "MonoisotopicMass", "AverageMass")
    
  } else {
    
    results <- tibble(UNIPROTKB = pull(tbl, accession_name))
    
  }

  results <- results %>% 
    left_join(dplyr::select(results_protein[[filelistnum]],
                            -c(Qvalue, 
                               filename,
                               monoiso_mass,
                               ave_mass,
                               fraction)), by = "UNIPROTKB")
  
  
  return(results)
  
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
  
  message("Getting GO subcellular locations...")
  
  temptbl <- tibble()
  
  temptbl <- tbl %>% 
    add_column(GO_term = NA) %>%
    add_column(GO_subcell_loc = NA)
  
  for (i in seq_along(tbl$`GO-ID`)) {
    
    if (is.na(temptbl$`GO-ID`[i]) == FALSE) {
    
    temptbl$GO_term[i] <- unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>% 
      trimws() %>% Term() %>% paste(collapse = "; ")
    
    temptbl$GO_subcell_loc[i] <- unlist(strsplit(temptbl$`GO-ID`[i], ";")) %>%
      trimws() %>% Term() %>% .[. %in% go_locs] %>% paste(collapse = "; ")
    
    }
  }
  
  return(temptbl)
}

getlocations <- function(resultslist) {
  
  # This function gets counts of membrane, cytosolic, and "both" proteoforms based on
  # GO terms pulled from UniProt for each unique accession number.
  
  counts <- tibble(filename = basename(names(resultslist)),
                   proteoform_count = NA,
                   cytosol_count = NA,
                   membrane_count = NA,
                   both_count = NA,
                   none_count = NA)
  
  for (i in seq_along(resultslist)) {
    
    # For every proteoform in each output, get the count of proteoforms whose GO terms
    # include "cytosol" OR "cytoplasm", "membrane", or BOTH. 
    # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!
    
    allAccession <- resultslist[[i]]$UNIPROTKB
    
    bothAccession <-  resultslist[[i]]$UNIPROTKB[str_detect(resultslist[[i]]$GO_subcell_loc,
                                                            c("cytosol|cytoplasm", "membrane"))]
    
    cytosolAccession <-  resultslist[[i]]$UNIPROTKB[str_detect(resultslist[[i]]$GO_subcell_loc,
                                                               c("cytosol|cytoplasm"))]
    
    membraneAccession <- resultslist[[i]]$UNIPROTKB[str_detect(resultslist[[i]]$GO_subcell_loc,
                                                               c("membrane"))]
    
    glue("{sum(cytosolAccession %in% bothAccession)} cytosolic proteoforms in 'both' for iteration {i}") %>% message
    
    glue("{sum(membraneAccession %in% bothAccession)} membrane proteoforms in 'both' for iteration {i}") %>% message
    
    # For cytosol_count and membrane_count, ONLY count the accessions which are NOT 
    # found in the list of accessions including both "cytosol|cytoplasm" and "membrane".
    # This prevents double-counting of proteoforms by localization
    
    counts$proteoform_count[i] <- resultslist[[i]] %>% .$UNIPROTKB %>% length()
    
    counts$cytosol_count[i] <- sum(!cytosolAccession %in% bothAccession)
    
    counts$membrane_count[i] <- sum(!membraneAccession %in% bothAccession)
    
    counts$both_count[i] <- length(bothAccession)
    
    counts$none_count[i] <- sum(!(allAccession %in% bothAccession) &
                                  !(allAccession %in% cytosolAccession) &
                                  !(allAccession %in% membraneAccession))
    
  }
  
  return(counts)
}

# Read Data Files -----------------------------------------------------------------------------

# Files are read from subdirectory called "input". Data files must
# be in csv or xlsx format and have a column of UniProt IDs whose 
# name includes the word "accession" somewhere. Case doesn't matter
# but spelling does.

setwd(filedir)

extension <- filelist %>% future_map(tools::file_ext)

if (length(unique(extension)) > 1) {
  
  stop("More than one kind of file. Try again.")
  
} else if (length(extension) == 0) {
  
  stop("No acceptable input files. Only tdReport, csv, and xlsx allowed.")
  
} else if (extension[[1]] == "csv") {
  
  message("Input is csv file, no proteoform data!")
  
} else if (extension[[1]] == "xlsx") {
  
  message("Input is xlsx file, no proteoform data!")

} else if (extension[[1]] == "tdReport") {
  
  message("Reading proteoform data from tdReport...")
  proteoformlist <- filelist %>% future_map2(., proteinlist, read_tdreport,
                                      fdr_cutoff = fdr,
                                      file_dir = filedir)
  tdreport_file <- TRUE
  
} else {
  print("What are you doing with your life?")
  stop()
}


# *** THE REST OF THE SCRIPT WILL ONLY EXECUTE IF THE INPUT IS 
# *** TDREPORT FORMAT!!

if (tdreport_file == TRUE) {
  
# **************************

names(proteoformlist) <- filelist

setwd(here())

# Access UniProt ------------------------------------------------------------------------------

message("Getting UniProt info from protein results...")

# keytypes <- UniProt.ws::keytypes(UPtaxon) %>% enframe
# columns <- UniProt.ws::columns(UPtaxon) %>% enframe

filelistnum <- as.list(seq_along(filelist))

results_proteoform <- proteoformlist %>% 
  future_map2(filelistnum,
              getuniprotinfo2,
              tdreport = tdreport_file)

if (tdreport_file == FALSE) {
  
  message("Adding masses according to UniProt sequences...")
  results_proteoform %<>% future_map(addmasses)
  
} else {
  
  message("Skipping addition of average and monoisotopic masses!")
  
}

results_proteoform[[length(results_proteoform)+1]] <- getlocations(results_proteoform)


# Output --------------------------------------------------------------------------------------

# An xlsx file with sheets corresponding to files in /data is written to
# the project directory

resultsname <- glue("{systime}_proteoform_results.xlsx")
resultsobjectname <- glue("{systime}_proteoform_results.rds")

names(results_proteoform) <- unlist(filelist)
names(results_proteoform)[length(results_proteoform)] <- "SUMMARY"

for (i in seq_along(names(results_proteoform))) {
  
  names(results_proteoform)[i] <- str_replace_all(names(results_proteoform[i]), "[:punct:]", "")
  names(results_proteoform)[i] %<>% stringr::str_trunc(28, "left") %>% paste(i, "_", ., sep = "")
  
}

setwd(here("output"))

results_proteoform %>%
  writexl::write_xlsx(path = resultsname)

results_proteoform %>%
  saveRDS(file = glue("rds/{resultsobjectname}"))

}

setwd(here())

