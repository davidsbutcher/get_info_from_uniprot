if(!require(here)) install.packages("here")

if(!require(svMisc)) install.packages("svMisc")

if(!require(furrr)) install.packages("furrr")

if(!require(Peptides)) install.packages("Peptides")

if(!require(magrittr)) install.packages("magrittr")

if(!require(writexl)) install.packages("writexl")

if(!require(readxl)) install.packages("readxl")

if(!require(readxl)) install.packages("openxlsx")

if(!require(RSQLite)) install.packages("RSQLite")

if(!require(tools)) install.packages("tools")

if(!require(progress)) install.packages("progress")

if(!require(tictoc)) install.packages("tictoc")

if(!require(glue)) install.packages("glue")

if(!require(RPushbullet)) install.packages("RPushbullet")

if(!require(ggpubr)) install.packages("ggpubr")

if(!require(tidyverse)) install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(GO.db)) BiocManager::install("GO.db")

if(!require(UniProt.ws)) BiocManager::install("UniProt.ws")

