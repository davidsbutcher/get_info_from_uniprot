
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

message("Retrieving info from UniProt...")

results <- tibble()

# Create a safe version of UniProt.ws which will not crash
# the whole damn program if it fails

accession_chunks <- chunk2(accession_numbers, 200)

safeselect <- safely(UniProt.ws::select)

glue("There are a total of {length(accession_chunks)} chunks") %>% 
  message

pb <- progress_bar$new(
  format = "  Accessing UniProt [:bar] :percent eta: :eta :spin",
  total = length(accession_chunks),
  clear = FALSE, width= 60)

for (i in seq_along(accession_chunks)) {
  
  pb$tick()
  
  suppressMessages(saferesults <- 
                     safeselect(UPtaxon,
                                keys = accession_chunks[[i]],
                                columns = c("ENTRY-NAME", "GENES",
                                            "PROTEIN-NAMES", "ORGANISM", 
                                            "ORGANISM-ID", "SEQUENCE", 
                                            "FUNCTION",
                                            "SUBCELLULAR-LOCATIONS", 
                                            "GO-ID"),
                                keytype = "UNIPROTKB"))
  
  if (is.null(saferesults[["result"]]) == TRUE) {
    
    glue("FAILURE! Could not access UniProt webservice for chunk {i}!") %>% 
      message
    
  } else {
    
    results <- union_all(results, saferesults[["result"]])
    
  }
  
}

UPdatabase <- results

saveRDS(results, 
        file = glue("input/{taxon_number}_full_UniProt_database.rds"))

rm(results)
        