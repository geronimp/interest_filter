#!/usr/bin/env Rscript
# Usage: basic_heatmap.R -f input_matrix.tsv -o output_heatmap.svg

library('pheatmap')
library("optparse")

################################################################################

option_list = list(
  	make_option(c("-f", "--input_matrix"), 
  				type="character", 
  				default=NULL, 
            	help="Input matrix file name", 
              	metavar="character"),
    make_option(c("-o", "--output_pdf"), 
    			type="character", 
    			default=NULL, 
              	help="Output svg name", 
              	metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$input_matrix)){
  print_help(opt_parser)
  stop("An input matrix must be specified", call.=FALSE)
}
if (is.null(opt$output_pdf)){
  print_help(opt_parser)
  stop("An output svg must be specified", call.=FALSE)
}

################################################################################
thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

PATH = dirname(thisFile())
ko_descriptions = paste(PATH, 'data/ko_descriptions.27-07-2017.tsv', sep='/')
data_matrix = read.delim(opt$input_matrix,
                         row.names=1)
ko2description     = read.delim(ko_descriptions,
                         quote="",
                         header = FALSE)

################################################################################

rownames(data_matrix) = paste(rownames(data_matrix), 
                              ko2description[match(rownames(data_matrix), ko2description$V1),]$V2,
                              sep="; ")

pheatmap(data_matrix, 
		 colorRampPalette(c('white', 'black'))(100),
		 border_color = 'gray50',
		 cellwidth    = 12,
		 cellheight   = 12,
		 fontsize     = 12,
     cluster_cols = FALSE,
     cluster_rows = FALSE,
     filename     = opt$output_pdf)
