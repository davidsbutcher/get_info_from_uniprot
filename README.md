
# Get UniProt Protein Info (GUPPI)

Process TDReports or flat files by getting information from the UniProt webservice and filtering by a selectable FDR value. Install the necessary packages, set input parameters and source "00_run_all.R" to process data.

## Input

All input parameters are given in 00_run_all.R, in the "Initialize Parameters" section. Acceptable input files are tdReport, xlsx, or csv format. Only one kind of file can be processed at a time. Data in xlsx or csv files must be a list of UniProt accession numbers with a column name that includes the word "Accession". Capitalization doesn't matter but spelling does. 

- `filedir <- c("Path to Folder or File")` Add the path of files to be processed. There are several options for adding input files: 
   
    - Enter the full path to the folder containing multiple files to be scanned, e.g. `filedir <- c("Z:/ICR/David Butcher/TDReports/")`. By default all subdirectories will be checked *unless* the directory has "deprecated" in its name. **Make sure to add the final forward slash** for directory names. 

    - Enter the full path to a single file, e.g. `filedir <- c("Z:/ICR/David Butcher/TDReports/Specific1.tdReport")`

    - Enter the full path to multiple single files, e.g. `filedir <- c("Z:/ICR/David Butcher/TDReports/Specific1.tdReport", "Z:/ICR/David Butcher/TDReports/Specific2.tdReport")`

- `fdr <- 0.01` This is the value used for the False Detection Rate cutoff. Defaults to 0.01 (1% FDR).

- `taxon_number <- 83333` UniProt taxon number for organism of interest. Defaults to 83333 for *E. coli* K12. Value for *Homo sapiens* is 9606. Taxons 83333 and 9606 are predownloaded (in /input), any other taxon number will take a minute to download.

- `go_locs_file <- "QuickGO_annotations_20190708.tsv"` This file will be used to determine which GO terms correspond to subcellular locations. This should not need to be changed.

- `plan(multisession(workers = 10))` This is a parameter needed by the furrr package. It designates the number of simultaenous R sessions used for data processing. Using 10 workers allows for 10 files to be processed simultaneously, but uses a lot of system resources. Drop this number on weaker systems or if <10 files are being processed. If processing fewer than 10 files, set the number of workers to the number of files to avoid system slowdown.

## Analysis

### Flat Files

The UniProt web service is queried for all UniProt accession numbers in the file using the package `UniProt.ws`. If a corresponding entry is found then protein name, organism, organism taxon ID, protein sequence, protein function, subcellular location, and any associated GO IDs are returned. Note that some of these values may not be found and come back as empty or NA. GO terms are obtained for all GO IDs using the `GO.db` package and terms corresponding to subcellular locations are saved in column "GO_subcellular_locations". Average and monoisotopic masses are determined using the `Peptides` package.

### TD Reports

A connection is established to the SQLite database in the TD Report using `RSQLite`. All protein isoform accession numbers are obtained, as well as the lowest Q value from among all hits for each  isoform and the name of the data file from which the lowest Q value hit was obtained. All isoforms with Q values which are missing or greater than the cutoff value are deleted.

The UniProt web service is queried for all remaining UniProt accession numbers using the package `UniProt.ws`. If a corresponding entry is found then protein name, organism, organism taxon ID, protein sequence, protein function, subcellular location, and any associated GO IDs are returned. Note that some of these values may not be found and come back as empty or NA. GO terms are obtained for all GO IDs using the `GO.db` package and terms corresponding to subcellular locations are saved in column "GO_subcellular_locations".  Average and monoisotopic masses are taken directly from the value in the TD Report.

For TD Report files, minimum Q value from among all hits, average and monoisotopic masses, and data file for lowest Q value hit are obtained for all proteoforms. Proteoforms with Q values or corresponding protein entry Q values above the FDR cutoff are deleted. UniProt info for every unique proteoform record number is copied from corresponding protein entries to avoid wasting time by querying the UniProt web service again.

## Output

Output files are saved to the output directory. Files are timestamped with the time the script was initialized or share the same name as the input file. 

* Protein results are saved to `output/YYYYMMDD_hhmmss_protein_results.xlsx` and `rds/YYYYMMDD_hhmmss_protein_results.rds`.
* Proteoform results are saved to `output/YYYYMMDD_hhmmss_proteoform_results.xlsx` and `rds/YYYYMMDD_hhmmss_protein_results.rds`.
  + Output xlsx files contain sheets corresponding to each input file specified in 00_run_all.R and a final summary sheet containing all input file names and counts of cytosolic/membrane proteins.
* Lists of all protein hits including UniProt accession number, Q value, data file, result type (i.e. tight absolute mass, find unexpected modifications, or biomarker) are saved to `output/allproteinhits` and share names with input files.
* Lists of all unique protein hits for each data file in a TDReport (intended for use in UpSet plots, pie charts, waffle plots, etc.) are saved to `output/proteinsbydatafile`.
* Mass histograms for identified proteins and proteoforms are saved to `/png` and `/pdf` in the corresponding formats.
* Workspace images (R objects) are saved to `output/workspace_images`.

### Pushbullet Notification

If your phone or web browser are set up for Pushbullet notification, a message can be sent to alert you of the script being finished. Check the function pbSetup in the RPushBullet package for more info.

