library(purrr)
library(vroom)
library(arrow)

outdir <- fs::dir_create("data")

move_files <- function(file) {
	file.rename(file, paste0("data/", fs::path_file(file)))
}

save_parquet <- function(file) {
  gtf_column_names = c("seqname", "source", "feature", 
		"start", "end", "score", "strand", "phase", "attribute") # standard column names for gtf files
  path <- fs::path_ext_remove(file) |> fs::path_ext_set("parquet") |> fs::path_file()
 if(grepl(file, pattern = "gtf"))
 	df <- vroom::vroom(file,comment="#", delim="\t", col_names = gtf_column_names)
 else
 	df <- vroom::vroom(file,comment="#", delim="\t")
  arrow::write_parquet(df,fs::path(outdir,path))
}

# WRITE OUTS ================================================================================
fs::dir_ls(outdir) |> fs::file_delete() # delete files present in the directory
fs::dir_ls("cache",regexp=".gz") |> walk(save_parquet) # convert .gz files to parquet
fs::dir_ls("cache", regexp=".h5ad") |> walk(move_files) # h5ad files will not be converted as this format is convenient for input into other programs directly
