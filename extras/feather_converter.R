rdsToConvert <- "input/XXXX.rds"

feathername <- 
  rdsToConvert %>% 
  basename() %>%
  fs::path_ext_remove() %>% 
  paste("input/", ., ".feather", sep = "")

rds <- readRDS(rdsToConvert)

# rds <- 
#   rds %>% as_tibble

feather::write_feather(rds, feathername)


TEST <- feather::read_feather("input/9606_full_UniProt_database.feather")
