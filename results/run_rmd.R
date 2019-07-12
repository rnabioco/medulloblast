#!/usr/bin/env Rscript

# run like this:
#Rscript --vanilla run_rmd.R file.Rmd n_cores

args = commandArgs(trailingOnly=TRUE)

rmarkdown::render(args[1], params = list(ncores = args[2]))
