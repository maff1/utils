#!/usr/bin/env Rscript

setUser = Sys.info()[["user"]]
LibDir = paste0("C:\\Users\\", setUser, "\\AppData\\Local\\R\\win-library\\4.2")
if(!dir.exists(LibDir)) {dir.create(LibDir)}

install.packages(setdiff(c("optparse", "dplyr", "openxlsx"), rownames(installed.packages())),
lib=LibDir,
repos="https://cran.ma.imperial.ac.uk",
dependencies=TRUE)

library("optparse")
library("dplyr")
library("openxlsx")

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="directory for files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="summary-stats.xlsx", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input dir).n", call.=FALSE)
}

f_format_mean_sd <- function(svar, sdig, ssep, sns) {
  paste(format(mean(svar, na.rm=T),digits=sdig, nsmall=sns), 
        format(sd(svar, na.rm=T),digits=sdig, nsmall=sns), sep = ssep)
}

f_format_min_max <- function(svar, sdig, ssep, sns) {
  paste(format(min(svar, na.rm=T),digits=sdig, nsmall=sns), 
        format(max(svar, na.rm=T),digits=sdig, nsmall=sns), sep = ssep)
}

ls_files <- list.files(path = opt$dir, pattern = ".csv", recursive = TRUE, full.names = F)
print(paste0(length(ls_files), " files found"))
print(basename(ls_files))
names(ls_files) <- substr(gsub("\\.\\/|\\/|\\.csv", "", ls_files), 1, 30)
  ls_df <- lapply(ls_files, function (x) { 
    df <- read.csv(x)
    out_df <- df %>%
     group_by(Source.File) %>%
     summarise(
       across(
         where(~ is.numeric(.x)),
           list(
       "min.max" = ~ f_format_min_max(svar=.x, sdig=2, sns=2, ssep=":"),
       "mean.sd" = ~ f_format_mean_sd(svar=.x, sdig=2, sns=2, ssep="+/-")
           ), 
         .names = "{.col}.{.fn}"),
       "N.total.peptides" = n(),
       "N.distinct.peptides" = n_distinct(Peptide)
     )
    out_df
    }
  )
openxlsx::write.xlsx(ls_df, asTable = T, file = opt$out)
print(paste0("summary-stats.xlsx", " in: ", getwd()))