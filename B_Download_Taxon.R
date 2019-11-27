
# Accessory script used to download info on all UniProt entries in a taxon.
# It is better to run this once and lookup info from the database then to
# access UniProt every time a TDReport is analyzed

# Packages ----------------------------------------------------------------

library(UniProt.ws)
library(purrr)
library(progress)
library(magrittr)
library(dplyr)
library(tibble)

# Functions ---------------------------------------------------------------

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# Get Accession Numbers from a Taxon ID -----------------------------------

accession_numbers <- UPtaxon@taxIdUniprots

# Download Data -----------------------------------------------------------

results <- tibble()

oldworkernum <- nbrOfWorkers()
plan(multisession(workers = 5))

# Create a safe version of UniProt.ws which will not crash
# the whole damn program if it fails

numberofchunks <- ceiling(length(accession_numbers)/100)

accession_chunks <- chunk2(accession_numbers, numberofchunks)

safeselect <- safely(UniProt.ws::select)

glue("There are {length(accession_chunks)} chunks") %>% 
  message

colsToQuery <- c("ENTRY-NAME", "GENES",
            "PROTEIN-NAMES", "ORGANISM", 
            "ORGANISM-ID", "SEQUENCE", 
            "FUNCTION",
            "SUBCELLULAR-LOCATIONS", 
            "GO-ID")

glue("Getting info from UniProt for taxon {taxon_number}..") %>%
  message

results_safe <- 
  future_map(accession_chunks, ~safeselect(UPtaxon,
                                           keys = .x,
                                           columns = colsToQuery,
                                           keytype = "UNIPROTKB"),
             .progress = TRUE
  )

results_safe %>%
  walk(~(if (is.null(.x[["result"]]) == TRUE) stop("NULL RESULT FOUND")) )

results <-
  results_safe %>%
  map(~(.x[["result"]])) %>%
  reduce(union_all) %>%
  as_tibble

# Output and Cleanup ------------------------------------------------------

UPdatabase <- results

saveRDS(results, 
        file = glue("input/{taxon_number}_full_UniProt_database.rds"))

feather::write_feather(results,
                       glue("input/{taxon_number}_full_UniProt_database.feather"))

rm(results)

plan(multisession(workers = oldworkernum %>% as.integer))
