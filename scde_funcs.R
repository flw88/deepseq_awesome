# Functions for 

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

ReadCellInfo <- function(filename){
  library(plyr)
  cell.info <- read.csv(filename, header=T)
  
  #Add some specific columns
  
  # broad_class_type
  new.vals <- c('gabaergic', 'gabaergic', 'gabaergic', 'distinct', 'nonneuronal', 'glutamatergic')
  names(new.vals)<- c("Gad2_Chrna6", "Gad2_Sepp1", "Gad2_Syt4", "Lars2_Kcnmb1", "Olig1", "Slc17a6")
  cell.info[['broad_class_type']] <- revalue(cell.info[['broad_class']], replace=new.vals)
  
  rownames(cell.info) <- cell.info[['rnaseq_profile_id']]
  
  return(cell.info)
}

# function for calculating ranking statistics from a DE results table
computeRankingStat <- function(result.table){
  if(!all(c('Z', 'p_value') %in% colnames(result.table))){
    stop("Function computeRankingStat expects matrix columns 'Z' and 'p_value'")
  }
  return(sign(result.table[,'Z'])*-1*log10(result.table[,'p_value']))
}