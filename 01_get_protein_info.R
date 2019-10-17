# Excel sheet names are limited to 30 characters.
# This script truncates the filenames after 30
# characters to use as Excel sheet names.

# Functions -----------------------------------------------------------------------------------

read_tdreport <- function(tdreport, fdr_cutoff = 0.01) {
  
  message(glue("\nEstablishing connection to {basename(tdreport)}..."))

  #Establish database connection. Keep trying until it works!
  #
   
  safe_dbConnect <- safely(dbConnect)
  
  safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)
  
  if (is.null(safecon[["result"]]) == TRUE) message("Connection failed, trying again!")
  
  iteration_num <- 1
  
  while (is.null(safecon[["result"]]) == TRUE) {
   
    iteration_num <- iteration_num + 1
    
    message(glue("\nTrying to establish database connection, attempt {iteration_num}"))
    safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:",
                              dbname = tdreport,
                              synchronous = NULL)
    
  }
  
  con <- safecon[["result"]]
  
  # Initialize the tibble where we'll be dropping our most 
  # important data, add NA columns for later use
  
  isoform_id <- con %>%
    RSQLite::dbGetQuery("SELECT Id, AccessionNumber 
                        FROM Isoform") %>%
    as_tibble %>% 
    add_column(Qvalue = NA) %>% 
    add_column(filename = NA)
  
  # Get Qvals and other info from "GlobalQualitativeConfidence"
  # table
  
  q_vals <- con %>%
    RSQLite::dbGetQuery("SELECT ExternalId, GlobalQvalue, HitId 
                        FROM GlobalQualitativeConfidence") %>%
    as_tibble %>% filter(.$ExternalId != 0)
  
  for (i in seq_along(isoform_id$Id)) {
    
    qValueExists <- any(!is.na(q_vals$GlobalQvalue[q_vals$ExternalId == isoform_id$Id[i]]))
    
    # If there is a Q value, get minimum Q value from among all Q values
    # for this Isoform ID then assign it to "Qvalue" column. Otherwise
    # assign value 1000000 so it is rejected when performing FDR cutoff
    
    if (qValueExists == TRUE) {
      
      min_q_val <-
        q_vals$GlobalQvalue[q_vals$ExternalId == isoform_id$Id[i]] %>%
        min(na.rm = TRUE)

      
    } else {
      
      print(glue("NO VALID Q VALUE FOUND, inserting garbage value in row {i}"))
      
      min_q_val <- 1000000
      
    }
    
    isoform_id$Qvalue[i] <- min_q_val
    
    # Get HitId from q_vals corresponding to lowest Q value hit
    # to be used to find information on data file name
    
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
      
      isoform_id$filename[i] <- con %>%
        RSQLite::dbGetQuery("SELECT Id, Name 
                          FROM DataFile") %>%
        as_tibble %>% 
        filter(Id == datafile_id) %>% 
        .$Name
      
    } else {
      
      isoform_id$filename[i] <- NA
      
    }
  }
  
  # Filter isoform_id by FDR cutoff value (default 0.01 or 1%)
  # and place results in "output" tibble
  
  output <- isoform_id %>% 
    filter(Qvalue < fdr_cutoff)
  
  setwd(here())
  
  # Close database connection and return output table
  
  dbDisconnect(con)
  
  message("read_tdreport Finished!")

  return(output)
  
}

read_tdreport_full <- function(tdreport, fdr_cutoff = 0.01) {
  
  # This function will return ALL hits above Q value threshold, not only the one
  # with the lowest Q value. 
  # Output is a tibble with all proteins hits below
  # FDR cutoff
  
  message(glue("\nEstablishing connection to {basename(tdreport)}..."))
  
  #Establish database connection. Keep trying until it works!
  
  safe_dbConnect <- safely(dbConnect)
  
  safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)
  
  if (is.null(safecon[["result"]]) == TRUE) message("Connection failed, trying again!")
  
  iteration_num <- 1
  
  while (is.null(safecon[["result"]]) == TRUE & iteration_num < 500) {
    
    iteration_num <- iteration_num + 1
    
    message(glue("\nTrying to establish database connection, attempt {iteration_num}"))
    safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:",
                              dbname = tdreport,
                              synchronous = NULL)
    
  }
  
  if (is.null(safecon[["result"]]) == TRUE) stop("Failed to connect using SQLite!")
  
  con <- safecon[["result"]]
  
  # Get relevant tables from TDReport using SQL
  
  {
    globalqualconf <- con %>%
      RSQLite::dbGetQuery("SELECT GlobalQvalue, HitId 
                        FROM GlobalQualitativeConfidence") %>%
      as_tibble
    
    datafile <- con %>%
      RSQLite::dbGetQuery("SELECT Id, Name 
                        FROM DataFile") %>%
      as_tibble %>% 
      dplyr::rename("DataFileId" = Id) %>% 
      dplyr::rename("filename" = Name)
    
    resultset <- con %>%
      RSQLite::dbGetQuery("SELECT Id, Name 
                        FROM ResultSet") %>%
      as_tibble %>% 
      dplyr::rename("ResultSetId" = Id) %>% 
      dplyr::rename("ResultSetName" = Name)
    
    isoform <- con %>%
      RSQLite::dbGetQuery("SELECT Id, AccessionNumber, Sequence 
                        FROM Isoform") %>%
      as_tibble %>% 
      dplyr::rename("IsoformId" = Id) %>% 
      dplyr::rename("SEQUENCE" = Sequence)
    
    
    bioproteoform <- con %>%
      RSQLite::dbGetQuery("SELECT IsoformId, ChemicalProteoformId
                        FROM BiologicalProteoform") %>%
      as_tibble
    
    hit <- con %>%
      RSQLite::dbGetQuery("SELECT Id, ResultSetId, DataFileId, ChemicalProteoformId
                        FROM Hit") %>%
      as_tibble %>% 
      dplyr::rename("HitId" = Id)
    
    
    entry <- con %>%
      RSQLite::dbGetQuery("SELECT Id, AccessionNumber
                        FROM Entry") %>%
      as_tibble %>% 
      dplyr::rename("EntryId" = Id)
    
    # Get Qvals and other info from "GlobalQualitativeConfidence"
    # table, put it into a new tibble. Remove all values for Q
    # values less than FDR cutoff
    
    q_vals <- con %>%
      RSQLite::dbGetQuery("SELECT Id, ExternalId, GlobalQvalue, HitId
                        FROM GlobalQualitativeConfidence") %>%
      as_tibble %>%
      filter(.$ExternalId != 0) %>%
      filter(.$GlobalQvalue <= fdr_cutoff) %>%
      dplyr::rename("EntryId" = ExternalId)
    
  }

  allproteinhits <- 
    left_join(hit, globalqualconf) %>% 
    filter(GlobalQvalue <= 0.01) %>% 
    left_join(datafile) %>%
    left_join(resultset) %>%
    left_join(bioproteoform) %>% 
    left_join(isoform) %>% 
    select(-c("HitId", "ChemicalProteoformId")) %>%
    select("AccessionNumber", "filename",
           "GlobalQvalue", "ResultSetName", everything()) %>% 
    drop_na()
  
  # allproteinhits <- 
  #   left_join(isoform, bioproteoform) %>%
  #   left_join(hit %>%
  #               filter(HitId %in% q_vals$HitId),
  #             by = "ChemicalProteoformId") %>%
  #   left_join(q_vals) %>% 
  #   dplyr::select(-c(ChemicalProteoformId, Id, EntryId, HitId)) %>% 
  #   left_join(datafile) %>%
  #   left_join(resultset) %>%
  #   drop_na
  
  output <- allproteinhits
  
  setwd(here())
  
  # Close database connection and return output table
  
  dbDisconnect(con)
  
  message("read_tdreport_full Finished!")
  
  return(output)
  
}

read_tdreport_byfilename <- function(tdreport, fdr_cutoff = 0.01) {
  
  # This function will return ALL hits above Q value threshold, not only the one
  # with the lowest Q value. 
  # Output is a tibble with all proteins split up by filename
  
  message(glue("\nEstablishing connection to {basename(tdreport)}..."))
  
  #Establish database connection. Keep trying until it works!
  
  safe_dbConnect <- safely(dbConnect)
  
  safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)
  
  if (is.null(safecon[["result"]]) == TRUE) message("Connection failed, trying again!")
  
  iteration_num <- 1
  
  while (is.null(safecon[["result"]]) == TRUE & iteration_num < 500) {
    
    iteration_num <- iteration_num + 1
    
    message(glue("\nTrying to establish database connection, attempt {iteration_num}"))
    safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:",
                              dbname = tdreport,
                              synchronous = NULL)
    
  }
  
  if (is.null(safecon[["result"]]) == TRUE) stop("Failed to connect using SQLite!")
  
  con <- safecon[["result"]]
  
  # Get relevant tables from TDReport using SQL
  {
    globalqualconf <- con %>%
      RSQLite::dbGetQuery("SELECT GlobalQvalue, HitId 
                        FROM GlobalQualitativeConfidence") %>%
      as_tibble
    
    datafile <- con %>%
      RSQLite::dbGetQuery("SELECT Id, Name 
                        FROM DataFile") %>%
      as_tibble %>% 
      dplyr::rename("DataFileId" = Id) %>% 
      dplyr::rename("filename" = Name)
    
    isoform <- con %>%
      RSQLite::dbGetQuery("SELECT Id, AccessionNumber 
                        FROM Isoform") %>%
      as_tibble %>% 
      dplyr::rename("IsoformId" = Id)
    
    bioproteoform <- con %>%
      RSQLite::dbGetQuery("SELECT IsoformId, ChemicalProteoformId
                        FROM BiologicalProteoform") %>%
      as_tibble
    
    hit <- con %>%
      RSQLite::dbGetQuery("SELECT Id, DataFileId, ChemicalProteoformId
                        FROM Hit") %>%
      as_tibble %>% 
      dplyr::rename("HitId" = Id)
    
    
    entry <- con %>%
      RSQLite::dbGetQuery("SELECT Id, AccessionNumber
                        FROM Entry") %>%
      as_tibble %>% 
      dplyr::rename("EntryId" = Id)
    
    # Get Qvals and other info from "GlobalQualitativeConfidence"
    # table, put it into a new tibble. Remove all values for Q
    # values less than FDR cutoff
    
    q_vals <- con %>%
      RSQLite::dbGetQuery("SELECT Id, ExternalId, GlobalQvalue, HitId
                        FROM GlobalQualitativeConfidence") %>%
      as_tibble %>%
      filter(.$ExternalId != 0) %>%
      filter(.$GlobalQvalue <= fdr_cutoff) %>%
      dplyr::rename("EntryId" = ExternalId)
    
  }
  
  allproteinhits <-
    left_join(isoform, bioproteoform) %>%
    left_join(hit %>%
                filter(HitId %in% q_vals$HitId),
              by = "ChemicalProteoformId") %>%
    left_join(q_vals) %>%
    dplyr::select(-c(ChemicalProteoformId, Id, EntryId, HitId)) %>%
    left_join(datafile) %>%
    drop_na
  
  proteinhitsbyfilename <- 
    allproteinhits %>% 
    dplyr::select(-c("DataFileId", "IsoformId", "GlobalQvalue")) %>% 
    pivot_wider(names_from = "filename",
                values_from = "AccessionNumber",
                values_fn = list(AccessionNumber = list))
  
  output <- proteinhitsbyfilename
  
  setwd(here())
  
  # Close database connection and return output table
  
  dbDisconnect(con)
  
  message("read_tdreport_byfilename Finished!")
  
  return(output)
  
}

read_tdreport_byfraction <- function(tdreport, fdr_cutoff = 0.01) {
  
  # This function will return ALL hits above Q value threshold, not only the one
  # with the lowest Q value. 
  # Output is a tibble with all proteins split up by filename
  
  message(glue("\nEstablishing connection to {basename(tdreport)}..."))
  
  #Establish database connection. Keep trying until it works!
  
  safe_dbConnect <- safely(dbConnect)
  
  safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:", dbname = tdreport)
  
  if (is.null(safecon[["result"]]) == TRUE) message("Connection failed, trying again!")
  
  iteration_num <- 1
  
  while (is.null(safecon[["result"]]) == TRUE & iteration_num < 500) {
    
    iteration_num <- iteration_num + 1
    
    message(glue("\nTrying to establish database connection, attempt {iteration_num}"))
    safecon <- safe_dbConnect(RSQLite::SQLite(), ":memory:",
                              dbname = tdreport,
                              synchronous = NULL)
    
  }
  
  if (is.null(safecon[["result"]]) == TRUE) stop("Failed to connect using SQLite!")
  
  con <- safecon[["result"]]
  
  # Get relevant tables from TDReport using SQL
  {
    datafile <- con %>%
      RSQLite::dbGetQuery("SELECT Id, Name 
                        FROM DataFile") %>%
      as_tibble %>% 
      dplyr::rename("DataFileId" = Id) %>% 
      dplyr::rename("filename" = Name)
    
    isoform <- con %>%
      RSQLite::dbGetQuery("SELECT Id, AccessionNumber 
                        FROM Isoform") %>%
      as_tibble %>% 
      dplyr::rename("IsoformId" = Id)
    
    bioproteoform <- con %>%
      RSQLite::dbGetQuery("SELECT IsoformId, ChemicalProteoformId
                        FROM BiologicalProteoform") %>%
      as_tibble
    
    hit <- con %>%
      RSQLite::dbGetQuery("SELECT Id, DataFileId, ChemicalProteoformId
                        FROM Hit") %>%
      as_tibble %>% 
      dplyr::rename("HitId" = Id)
    
    
    entry <- con %>%
      RSQLite::dbGetQuery("SELECT Id, AccessionNumber
                        FROM Entry") %>%
      as_tibble %>% 
      dplyr::rename("EntryId" = Id)
    
    # Get Qvals and other info from "GlobalQualitativeConfidence"
    # table, put it into a new tibble. Remove all values for Q
    # values less than FDR cutoff
    
    q_vals <- con %>%
      RSQLite::dbGetQuery("SELECT Id, ExternalId, GlobalQvalue, HitId
                        FROM GlobalQualitativeConfidence") %>%
      as_tibble %>%
      filter(.$ExternalId != 0) %>%
      filter(.$GlobalQvalue <= fdr_cutoff) %>%
      dplyr::rename("EntryId" = ExternalId)
    
  }
  
  allproteinhits <- 
    left_join(isoform, bioproteoform) %>%
    left_join(hit %>%
                filter(HitId %in% q_vals$HitId),
              by = "ChemicalProteoformId") %>%
    left_join(q_vals) %>% 
    dplyr::select(-c(ChemicalProteoformId, Id, EntryId, HitId)) %>% 
    left_join(datafile) %>% 
    drop_na %>%
    addfraction()
  
  proteinhitsbyfraction <- 
    allproteinhits %>% 
    dplyr::select(-c("DataFileId", "IsoformId", "GlobalQvalue", "filename")) %>% 
    pivot_wider(names_from = "fraction",
                values_from = "AccessionNumber",
                values_fn = list(AccessionNumber = list))
  
  output <- proteinhitsbyfraction
  
  setwd(here())
  
  # Close database connection and return output table
  
  dbDisconnect(con)
  
  message("read_tdreport_byfilename Finished!")
  
  return(output)
  
}

getuniprotinfo <- function(tbl, taxon = NULL, tdreport = TRUE) {
  
  message("Retrieving info from UniProt...")
  
  # Find column in the input tibble which has "accession"
  # in it and use it to get info from UniProt
  
  accession_name <- grep("accession", names(tbl),
                         ignore.case = TRUE, value = TRUE)
  
  # If the data is coming from a TDreport, add a filename column
  # and Qvalue column
  
  if (tdreport == TRUE) {
    
    results <- dplyr::select(tbl, -c(1))
    names(results) <- c("UNIPROTKB", "Qvalue", "filename")
    
  } else {
    results <- tibble(UNIPROTKB = pull(tbl, accession_name))
  }
  
  results_temp <- tibble()
  
  # Prepare progress bar
  
  pb <- progress_bar$new(
    format = "Accessing UniProt [:bar] :percent eta: :eta :spin",
    total = length(pull(tbl, accession_name)),
    clear = FALSE, width= 60)
  
  # Create a safe version of UniProt.ws which will not crash
  # the whole damn program if it fails
  
  safeselect <- safely(UniProt.ws::select)
  
  for (i in seq_along(pull(tbl, accession_name))) {
    
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
      union_all(results_temp, results_newline)  
    
    pb$tick()
    
  }
  
  results_after_join <- left_join(results, results_temp)
  return(results_after_join)
  
  message("Finished retrieving info from UniProt!")
  
}

getuniprotinfo2 <- function(listofproteins, taxon = NULL, database = NULL, tdrep = TRUE) {
  
  message("Retrieving info from UniProt...")
  
  # Find column in the input tibble which has "accession"
  # in it and use it to get info from UniProt
  
  accession_name <- grep("accession", names(listofproteins),
                         ignore.case = TRUE, value = TRUE)
  
  if (tdrep == TRUE) {
    
    listofproteins <- listofproteins %>% 
      dplyr::select(-c(Id), "UNIPROTKB" := !!accession_name)
    
  } else {
    
    listofproteins <- tibble(UNIPROTKB = pull(listofproteins, accession_name))
    
  }
  
  left_join(listofproteins, database)
  
}


addmasses <- function(tbl) {
  
# This function uses the Peptides package to determine
# average and monoisotopic masses based on protein sequence.
# These values diverge from those in the TDreport by <0.00005 Da.

  mutate(tbl, monoiso_mass = mw(tbl$SEQUENCE, monoisotopic = TRUE),
         ave_mass = mw(tbl$SEQUENCE, monoisotopic = FALSE))
  
}

addfraction <- function(tbl) {
  
# This function attempts to parse the filenames to extract information
# about the fraction that a raw file corresponds to. This is only useful
# for GELFrEE/PEPPI/other fractionated data
  
  tbl %>% 
    mutate(
      fraction = case_when(
        stringr::str_detect(filename, "(?i)(?<=gf|peppi|frac|fraction|f|f_)[0-9]{1,2}") == TRUE ~
          stringr::str_extract(filename,
                               "(?i)(?<=gf|peppi|frac|fraction|f|f_)[0-9]{1,2}"),
        stringr::str_detect(filename, "(?i)(?<=gf|peppi|frac|fraction|f|f_)[0-9]{1,2}") == FALSE ~ "NA"
      )
    )
  
}

getGOterms <- function(tbl, file_list) {

  message(glue("Getting GO subcellular locations for {basename(file_list)}..."))
  
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
  
  # This function gets counts of membrane, cytosolic, and "both" proteins based on
  # GO terms pulled from UniProt for each unique accession number.
  
  counts <- tibble(filename = basename(names(resultslist)),
                   protein_count = NA,
                   cytosol_count = NA,
                   membrane_count = NA,
                   both_count = NA,
                   none_count = NA)
  
  for (i in seq_along(resultslist)) {

    # For every protein in each output, get the count of proteins whose GO terms
    # include "cytosol" OR "cytoplasm", "membrane", or BOTH. 
    # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!
    
    allAccession <- resultslist[[i]]$UNIPROTKB
    
    bothAccession <-  resultslist[[i]]$UNIPROTKB[str_detect(resultslist[[i]]$GO_subcell_loc,
                                                            c("cytosol|cytoplasm", "membrane"))]
    
    cytosolAccession <-  resultslist[[i]]$UNIPROTKB[str_detect(resultslist[[i]]$GO_subcell_loc,
                                                               c("cytosol|cytoplasm"))]
    
    membraneAccession <- resultslist[[i]]$UNIPROTKB[str_detect(resultslist[[i]]$GO_subcell_loc,
                                                               c("membrane"))]
    
    glue("{sum(cytosolAccession %in% bothAccession)} cytosolic proteins in 'both' for iteration {i}") %>%
      message
    
    glue("{sum(membraneAccession %in% bothAccession)} membrane proteins in 'both' for iteration {i}") %>%
      message

    # For cytosol_count and membrane_count, ONLY count the accessions which are NOT 
    # found in the list of accessions including both "cytosol|cytoplasm" and "membrane".
    # This prevents double-counting of proteins by localization
    
    counts$protein_count[i] <- resultslist[[i]] %>% .$UNIPROTKB %>% length()
    
    counts$cytosol_count[i] <- sum(!cytosolAccession %in% bothAccession)
    
    counts$membrane_count[i] <- sum(!membraneAccession %in% bothAccession)
    
    counts$both_count[i] <- length(bothAccession)
    
    counts$none_count[i] <- sum(!(allAccession %in% bothAccession) &
                                  !(allAccession %in% cytosolAccession) &
                                  !(allAccession %in% membraneAccession))
  
  }
  
  return(counts)
}

getlocations_protein <- function(resultslist) {
  
  # This function gets counts of membrane, cytosolic, and "both" proteoforms based on
  # GO terms pulled from UniProt for each unique accession number.
  
  counts <- tibble(filename = basename(names(resultslist)),
                   protein_count = NA,
                   cytosol_count = NA,
                   membrane_count = NA,
                   periplasm_count = NA,
                   NOTA_count = NA)
  
  for (i in seq_along(resultslist)) {
    
    tempresults <- tibble()
    
    # For every proteoform in each output, get the count of proteoforms whose GO terms
    # include "cytosol" OR "cytoplasm", "membrane", or BOTH. 
    # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!
    tempresults <- 
      resultslist[[i]] %>% 
      mutate(cytosol = str_detect(resultslist[[i]]$GO_subcell_loc,
                                  c("cytosol|cytoplasm|ribosome"))) %>%
      mutate(membrane = str_detect(resultslist[[i]]$GO_subcell_loc, 
                                   c("membrane")) &
               !str_detect(resultslist[[i]]$GO_subcell_loc, 
                           c("membrane-bounded periplasmic space"))) %>% 
      mutate(periplasm = str_detect(resultslist[[i]]$GO_subcell_loc,
                                    c("periplasm")))
    
    # For cytosol_count and membrane_count, ONLY count the accessions which are NOT 
    # found in the list of accessions including both "cytosol|cytoplasm" and "membrane".
    # This prevents double-counting of proteoforms by localization
    
    counts$protein_count[i] <- tempresults %>% .$UNIPROTKB %>% length()
    
    counts$cytosol_count[i] <- sum(tempresults$cytosol)
    
    counts$membrane_count[i] <- sum(tempresults$membrane)
    
    counts$periplasm_count[i] <- sum(tempresults$periplasm)
    
    counts$NOTA_count[i] <- tally(tempresults %>% 
                                    filter(cytosol == FALSE) %>% 
                                    filter(membrane == FALSE) %>% 
                                    filter(periplasm == FALSE)) %>%
      .$n
    
  }
  
  return(counts)
}

getlocations_bydatafile <- function(resultslist) {
  
  # This function gets counts of membrane, cytosolic, and "both" proteoforms based on
  # GO terms pulled from UniProt for each unique accession number.
  
  counts <- tibble()
  
  for (i in seq_along(resultslist)) {
    
    tempresults <- tibble()
    
    # For every proteoform in each output, get the count of proteoforms whose GO terms
    # include "cytosol" OR "cytoplasm", "membrane", or BOTH. 
    # WE DO NOT DIFFERENTIATE BETWEEN MEMBRANE TYPES!
    tempresults <- 
      resultslist[[i]] %>% 
      mutate(tdreport_name = names(resultslist)[[i]]) %>% 
      mutate(cytosol = str_detect(resultslist[[i]]$GO_subcell_loc,
                                  c("cytosol|cytoplasm|ribosome"))) %>%
      mutate(membrane = str_detect(resultslist[[i]]$GO_subcell_loc, 
                                   c("membrane")) &
               !str_detect(resultslist[[i]]$GO_subcell_loc, 
                           c("membrane-bounded periplasmic space"))) %>% 
      mutate(periplasm = str_detect(resultslist[[i]]$GO_subcell_loc,
                                    c("periplasm"))) %>% 
      mutate(NOTA = !str_detect(resultslist[[i]]$GO_subcell_loc,
                                c("cytosol|cytoplasm|ribosome")) &
               !str_detect(resultslist[[i]]$GO_subcell_loc, 
                           c("membrane")) &
               !str_detect(resultslist[[i]]$GO_subcell_loc, 
                           c("membrane-bounded periplasmic space")) &
               !str_detect(resultslist[[i]]$GO_subcell_loc,
                           c("periplasm"))
      )
    
    # 
    
    tempresultssummary <- 
      tempresults %>%
      group_by(tdreport_name, filename, fraction) %>%
      summarize(
        protein_count = n(),
        cytosol_count = sum(cytosol),
        membrane_count = sum(membrane),
        periplasm_count = sum(periplasm),
        NOTA_count = sum(NOTA)
      )
    
    counts <- 
      union_all(tempresultssummary, counts)
    
  }
  
  return(counts)
}

savePLBF <- function(input_tbbl, filename) {
  
  # Save proteinlist by filename
  
  workbook <- createWorkbook()
  
  shortfilename <- 
    filename %>%
    basename %>%
    file_path_sans_ext %>% 
    str_trunc(30, "left")
  
  workbook %>% addWorksheet(shortfilename)
  
  for (i in seq_along(input_tbbl)) {
    
    workbook %>% writeData(shortfilename,
                           names(input_tbbl)[[i]],
                           startCol = i,
                           startRow = 1)
    
    
    workbook %>% writeData(shortfilename,
                           pull(input_tbbl, i) %>% as_vector %>% unique,
                           startCol = i,
                           startRow = 2)
    
  }
  
  return(workbook)
  
}

savePLBFrac <- function(input_tbbl, filename) {
  
  # Save proteinlist by fraction
  
  workbook <- createWorkbook()
  
  shortfilename <- 
    filename %>%
    basename %>%
    file_path_sans_ext %>% 
    str_trunc(30, "left")
  
  workbook %>% addWorksheet(shortfilename)
  
  for (i in seq_along(input_tbbl)) {
    
    workbook %>% writeData(shortfilename,
                           names(input_tbbl)[[i]],
                           startCol = i,
                           startRow = 1)
    
    
    workbook %>% writeData(shortfilename,
                           pull(input_tbbl, i) %>% as_vector %>% unique,
                           startCol = i,
                           startRow = 2)
    
  }
  
  return(workbook)
  
}

# Read Data Files -----------------------------------------------------------------------------

setwd(filedir)

extension <- filelist %>% map(tools::file_ext)

if (length(unique(extension)) > 1) {
  
  stop("More than one kind of input file. Try again.")
  
} else if (length(extension) == 0) {
  
  stop("No acceptable input files. Only tdReport, csv, or xlsx are allowed.")
  
} else if (extension[[1]] == "csv") {
  
  message("Reading csv files...")
  proteinlist  <- filelist %>% map(read_csv)
  tdreport_file <- FALSE
  
} else if (extension[[1]] == "xlsx") {
  
  message("Reading xlsx files...")
  proteinlist  <- filelist %>% map(read_xlsx)
  tdreport_file <- FALSE
  
} else if (extension[[1]] == "tdReport") {
  
  message("Reading protein data from tdReport...")
  proteinlist <- filelist %>% future_map(read_tdreport, fdr_cutoff = fdr)
  
  message("Reading full protein data from tdReport...")
  proteinlistfull <- filelist %>% 
    future_map(read_tdreport_full, fdr_cutoff = fdr)
  
  message("Reading protein data by file name from tdReport...")
  proteinlistbyfilename <- filelist %>% 
    future_map(read_tdreport_byfilename, fdr_cutoff = fdr)
  
  message("Attempting to read protein data by fraction from tdReport...")
  proteinlistbyfraction <- filelist %>% 
    future_map(read_tdreport_byfraction, fdr_cutoff = fdr)

  tdreport_file <- TRUE
  
} else {
  print("What are you doing with your life?")
  stop()
}

names(proteinlist) <- filelist

setwd(here())

# Access UniProt ------------------------------------------------------------------------------

if (glue("input/{taxon_number}_full_UniProt_database.rds") %>% 
    file.exists == TRUE) {
  
  glue("Found UniProt database for taxon {taxon_number}, loading") %>% 
    message
  
  UPdatabase <- readRDS(glue("input/{taxon_number}_full_UniProt_database.rds"))
  
} else {
  
  glue("DID NOT FIND UniProt database for taxon {taxon_number}, DOWNLOADING...") %>% 
    message
  
   source("B_Download_Taxon.R")
  
}

results_protein <- proteinlist %>%
  future_map(getuniprotinfo2,
             taxon = taxon_number,
             database = UPdatabase,
             tdrep = tdreport_file) %>% 
  future_map2(filelist, getGOterms, .progress = TRUE) %>% 
  future_map(addmasses) %>% 
  future_map(addfraction)

proteinlistfull <- 
  proteinlistfull %>%
  future_map(addmasses) %>%
  future_map(addfraction)
  
results_protein[[length(results_protein)+1]] <- getlocations_protein(results_protein)

names(results_protein) <- unlist(filelist) %>% basename()

## This chunk will get protein location counts by datafile (.raw file) instead
## of for each TDreport in the analysis

protein_counts_bydatafile <- 
  results_protein[1 : length(results_protein) - 1] %>% 
  getlocations_bydatafile()

# Output --------------------------------------------------------------------------------------

# An xlsx file with sheets corresponding to files in /data is written to
# the project directory

# If systime doesn't exist (i.e. script is not being run from 00_run_all)
# create the systime variable

if (exists("systime") == FALSE) systime <- format(Sys.time(), "%Y%m%d_%H%M%S")

resultsname <- glue("{systime}_protein_results.xlsx")
resultsobjectname <- glue("{systime}_protein_results.rds")

names(results_protein) <- unlist(filelist) %>% basename()
names(results_protein)[length(results_protein)] <- "SUMMARY"

for (i in seq_along(names(results_protein))) {

  names(results_protein)[i] <- 
    str_replace_all(names(results_protein[i]), "[:punct:]", "")
  
  names(results_protein)[i] %<>% 
    stringr::str_trunc(28, "left") %>% paste(i, "_", ., sep = "")

}

if (dir.exists("output") == FALSE) dir.create("output")

setwd(here("output"))


# Save protein results ----------------------------------------------------

results_protein %>%
  writexl::write_xlsx(path = resultsname)

if (tdreport_file == TRUE) {
  
  # Save protein results, all hits ------------------------------------------
  
  names(proteinlistfull) <- unlist(filelist) %>% basename
  
  if (dir.exists("allproteinhits") == FALSE) dir.create("allproteinhits")
  
  PLFnameslist <- 
    filelist %>%
    map(basename) %>%
    map(file_path_sans_ext) %>%
    glue_data("allproteinhits/{.}.xlsx") %>% 
    as.list
  
  map2(proteinlistfull, PLFnameslist, ~write_xlsx(.x, path = .y))
  
  # Save list of proteins by filename ---------------------------------------
  
  if (dir.exists("proteinsbydatafile") == FALSE) dir.create("proteinsbydatafile")
  
  PLBFlist <- 
    map2(proteinlistbyfilename, filelist, savePLBF)
  
  PLBFnameslist <- 
    filelist %>%
    map(basename) %>%
    map(file_path_sans_ext) %>%
    glue_data("proteinsbydatafile/{.}.xlsx") %>% 
    as.list
  
  map2(PLBFlist, PLBFnameslist, ~saveWorkbook(.x, .y, overwrite = TRUE))
  
  # Save list of proteins by fraction ---------------------------------------
  
  if (dir.exists("proteinsbyfraction") == FALSE) dir.create("proteinsbyfraction")
  
  PLBFrac_list <- 
    map2(proteinlistbyfraction, filelist, savePLBFrac)
  
  PLBFrac_nameslist <- 
    filelist %>%
    map(basename) %>%
    map(file_path_sans_ext) %>%
    glue_data("proteinsbyfraction/{.}.xlsx") %>% 
    as.list
  
  map2(PLBFrac_list, PLBFrac_nameslist, ~saveWorkbook(.x, .y, overwrite = TRUE))
  

  # Save counts by datafile -------------------------------------------------

  if (dir.exists("countsbydatafile") == FALSE) dir.create("countsbydatafile")
  
  protein_counts_bydatafile %>% 
    writexl::write_xlsx(path = glue("countsbydatafile/{systime}_counts_bydatafile.xlsx"))
  
}

# Save RDS ----------------------------------------------------------------

if (dir.exists("rds") == FALSE) dir.create("rds")

results_protein %>%
  saveRDS(file = glue("rds/{resultsobjectname}"))

setwd(here())
