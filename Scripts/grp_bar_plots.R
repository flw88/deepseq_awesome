#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

#Make some barplots of the cell annotations so that we know what we have at a glance

#Hack for collaboration. Given a list of possible paths, pick the first one which exists
detectDirs <- function(dir.list){
  for(cd in dir.list){
    if(file.exists(cd)){
      ret.dir <- cd
      break;
    }
  }
  return(ret.dir)
}

# function for parsing a filter string in a contrast spec file
ParseFilterStr <- function(filter.str){
  out <- list()
  toks <- strsplit(filter.str,";", fixed=TRUE)[[1]]
  for(tok in toks){
    tt <- strsplit(tok, "=", fixed=TRUE)[[1]]
    key <- tt[[1]]
    val <- strsplit(tt[[2]], ",", fixed=TRUE)[[1]]
    out[[key]] <- val
  }
  return(out)
}

# function for selecting samples using a sample.info dataframe and the output of parse.filter.str()
SelectSamples <- function(sample.info, fixed.cols, verbose=F){
  bool.use.samples <- rep(TRUE, nrow(sample.info))
  for(fc in names(fixed.cols)){
    fixed.vals <- fixed.cols[[fc]]
    if(verbose){
      print(sprintf('%s == %s', fc, fixed.vals))
    }
    bool.use.samples <- bool.use.samples & sample.info[,fc] %in% fixed.vals
  }
  return(bool.use.samples)
}

MakeDateStr <- function(){
  return(format(Sys.time(), "%Y%m%d.%H%M%S"))
}

library(ggplot2)


print(sprintf("Started grp_hist_plots.R at %s", Sys.time()))

home.dirs <- c("/home/local/users/flwu/teamawesome/", "/results/G4017_Users/team_awesome/",
               "/Users/Felix/GitHub/")
home.dir <- detectDirs(home.dirs)

data.dir <- file.path(home.dir, "datasets")
scrna.dir <- file.path(data.dir, "scRNASeq")
run.dir <- file.path(home.dir, "deepseq_awesome")

cell.info <- read.csv(file.path(scrna.dir, "columns-cells.csv"), header=T)
rownames(cell.info) <- cell.info[['rnaseq_profile_id']]

targ.cols <- c('subclass', 'broad_class', 'sampling_region', 'slice_index')
stopifnot(all(targ.cols %in% names(cell.info)))

for(cur.col in targ.cols){
  cell.info[[cur.col]] <- as.factor(cell.info[[cur.col]])
  cur.fn <- file.path(home.dir, 'results', paste(cur.col, 'barplot', 'pdf', sep='.'))
  cur.plot <- ggplot(data=cell.info, mapping=aes_string(cur.col)) + geom_bar() +
    coord_flip()
  ggsave(filename = cur.fn, plot = cur.plot, width=11, height=9)
}