#!/usr/bin/env Rscript

# copied from https://github.com/Cibiv/trumicount/blob/master/install.R

doc <- '
Usage: install [ --install-program INSTALL ] [ --prefix PREFIX ]
Options:
--install-program INSTALL    `install` program to use [Default: install]
--prefix PREFIX              Installation prefix [Default: /usr]
'
ARGS <- docopt::docopt(doc, strip_names=TRUE)

system(paste(ARGS$`install-program`, "minionqc", paste0(ARGS$prefix, "/bin/"), sep=" "))