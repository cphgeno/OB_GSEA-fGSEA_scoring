suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(pheatmap)
  library(fgsea)
  library(qusage)
})

# source helper functions
get_script_dir <- function() {
  commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value, into = c("key", "value"), sep = "=", fill = "right"
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value) %>%
    dirname(.)
}
source(file.path(get_script_dir(), "fgsea_helper_functions.R"))

# parse options
option_list <- list(
  make_option("--output_dir", "-o",
    dest="output_dir", type="character",
    help="output directory where files will be saved"),
  make_option("--name",
    type = "character",
    help = "name of the dataset"),
  make_option("--preprocessing.sampleranks",
    dest="ranks_df_sample", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.samplediffranks",
    dest="ranks_df_diff", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.referenceranks",
    dest="reference_df", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.meta",
    dest="metadata", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--preprocessing.genesets",
    dest="genesets", type="character",
    help="input directory where datasets and genesets are found"),
  make_option("--input_type", # either RankExpr or DeltaCentroid
    dest="input_type", type="character",
    help="input_type depending on ranked input")
)
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


# --- Required argument check ---
defined_opts <- sapply(option_list, function(x) x@dest)
missing_opts <- defined_opts[vapply(opts[defined_opts], is.null, logical(1))]
if (length(missing_opts) > 0) {
  print_help(opt_parser)
  stop(paste("Missing required option(s):",
             paste(missing_opts, collapse=", ")),
       call.=FALSE)
}

TOOL_COLOURS <- list(
  fgsea_DeltaCentroid = "#044be6",
  fgsea_RankExpr = "#AEC7E7"
)

# check analysis type to fetch correct input
if (opts$input_type == "RankExpr"){
  ranks_df <- opts$ranks_df_sample
  tool_colour = TOOL_COLOURS[[paste0("fgsea_", opts$input_type)]]
} else {
  ranks_df <- opts$ranks_df_diff
  tool_colour = TOOL_COLOURS[[paste0("fgsea_", opts$input_type)]]
}

input_df <- read.table(ranks_df, sep = '\t', header = TRUE, check.names = FALSE)
input_df <- column_to_rownames(input_df, var = 'Geneid')
metadata_df <- read.table(opts$metadata, sep = '\t', header = TRUE, check.names = FALSE, fill = TRUE, comment.char = "")
metadata_df <- metadata_df[rowSums(is.na(metadata_df)) < ncol(metadata_df), ]
genesets_list <- read.gmt(opts$genesets)

fgsea_wrapper(opts$output_dir, input_df, metadata_df, genesets_list, opts$name, opts$input_type, tool_colour)
