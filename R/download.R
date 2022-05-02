# This file downloads all data files from GTEx
library(fs) 
library(XML)
library(purrr)

# strsplit_grepl and get_long_read_data_files are used to parse out the GTEx file names.
strsplit_grepl <- function(x, pattern, extension) { 
	temp = unlist(strsplit(x, pattern))
	temp[unlist(lapply(temp, FUN=grepl, pattern = extension))]
}

get_file_names <- function(files, pattern1, pattern2="</Key>", extension) {
	strsplit_grepl(files, pattern1, extension) |> strsplit_grepl(pattern2, extension)
}

options(timeout=1800) # download timeout

fs::dir_create("cache") # Create directory for downloads/intermediate files
download.file("https://storage.googleapis.com/gtex_analysis_v9/", dest = "cache/files.xml") # download intermediate file with file names into cache
file_names_long_read = get_file_names(readLines("cache/files.xml"), pattern1 = "long_read_data/", extension = ".gz") # extract file names of long read ata
file_names_short_read = get_file_names(readLines("cache/files.xml"), pattern1 = "snrna_seq_data/", extension = ".h5ad") # extract file names of long read ata

URLs = c(paste0("https://storage.googleapis.com/gtex_analysis_v9/long_read_data/", 
	file_names_long_read), paste0("https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/", file_names_short_read))
path = paste0("cache/", c(file_names_long_read, file_names_short_read)) # Create output file names
walk2(URLs, path, download.file) # Download files 