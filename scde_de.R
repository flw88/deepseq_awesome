#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

#Perform differential expression analysis, using SCDE

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

library(scde)
library(monocle)

print(sprintf("Started scde_de.R at %s", Sys.time()))

home.dirs <- c("/home/local/users/flwu/teamawesome/", "/results/G4017_Users/team_awesome/")
home.dir <- detectDirs(home.dirs)

data.dir <- file.path(home.dir, "datasets")
scrna.dir <- file.path(data.dir, "scRNASeq")
run.dir <- file.path(home.dir, "deepseq_awesome")

count.df <- read.csv(file.path(scrna.dir, "count_table.ceiling.csv"), header=T, row.names=1)
gene.info <- read.csv(file.path(scrna.dir, "rows-genes.csv"), row.names=1, header=T, quote='"')
cell.info <- read.csv(file.path(scrna.dir, "columns-cells.csv"), header=T)
rownames(cell.info) <- cell.info[['rnaseq_profile_id']]

spec.table <- read.delim(file.path(run.dir, "scde_de.specs.txt"), header=T, sep='\t', comment.char='#', as.is=T, colClasses='character')

n.cores <- 16

#Change id names
names(count.df) <- gsub("^X", replacement = '', names(count.df))

#Order info dfs
stopifnot(all(rownames(count.df) %in% rownames(gene.info)))
gene.info <- gene.info[rownames(count.df),]
stopifnot(all(names(count.df) %in% rownames(cell.info)))
cell.info <- cell.info[names(count.df),]

all.de.results <- list()
for(i in 1:nrow(spec.table)){
  #Read in, parse spec table info
  spec.row <- spec.row <- spec.table[i,]
  group.col <- spec.row[['group.col']]
  class.pair <- c(spec.row[['case']], spec.row[['control']])
  
  fixed.cols <- ParseFilterStr(spec.row[['filters']])
  fixed.cols[[group.col]] <- class.pair # Ensure there aren't any NAs
  keep.cells <- SelectSamples(cell.info, fixed.cols)
  
  if((spec.row[['exclude']] != "") && !is.null(spec.row[['exclude']])) {
    exclude.cols <- parse.filter.str(spec.row[['exclude']]) # Exclude samples in 'exclude'
    exclude.cells <- SelectSamples(cell.info, exclude.cols)
    keep.cells <- keep.cells & !(exclude.cells) # Apply 'exclude' 
  }
  
  #Get only the indicated samples/cells
  cur.cell.info <- cell.info[keep.cells,]
  cur.count.m <- as.matrix(count.df[,keep.cells])
  cur.count.m <- apply(cur.count.m, 2, function(x) {storage.mode(x) <- 'integer'; x}) # Recast as integer
  
  batch.col <- spec.row[['batch.col']]  # which variables do we want to use as batch vars?
  if(batch.col == '' || batch.col == '.'){
    batch.factor <- NULL
  } else {
    batch.factor <- factor(cur.cell.info[[batch.col]], levels=unique(cur.cell.info[[batch.col]]))
    names(batch.factor) <- rownames(cur.cell.info)
  }
  
  out.base <- spec.row[['out.base']]
  
  print(sprintf("Started %s comparison at %s", out.base, Sys.time()))
  
  # Get factor of DE group
  group.factor <- factor(cur.cell.info[[group.col]], levels=class.pair)
  names(group.factor) <- rownames(cur.cell.info)
  
  # Clean count frame
  cur.clean.m <- clean.counts(cur.count.m, min.lib.size=1000, min.reads = 1, min.detected = 1)
  
  # Calculate error model
  err.mod <- scde.error.models(counts = cur.clean.m, groups = group.factor, n.cores = n.cores,
                               threshold.segmentation = TRUE, save.crossfit.plots = TRUE, save.model.plots = TRUE, verbose = 1)
  
  # filter out cells that don't show positive correlation with the expected expression magnitudes (very poor fits)
  valid.cells <- err.mod$corr.a > 0
  print(sprintf("Keeping %s of %s total cells.", sum(valid.cells), length(valid.cells)))
  err.mod <- err.mod[valid.cells,]
  group.factor <- group.factor[valid.cells]
  
  # estimate gene expression prior
  e.prior <- scde.expression.prior(models = err.mod, counts = cur.clean.m, length.out = 400, show.plot = FALSE)
  if(is.null(batch.factor)){
    ediff <- scde.expression.difference(err.mod, count.df, e.prior, groups = group.factor, n.randomizations = 100,
                                        n.cores = n.cores, verbose = 1)
  } else {
    batch.factor <- batch.factor[valid.cells]
    ediff <- scde.expression.difference(err.mod, count.df, e.prior, groups = group.factor, n.randomizations = 100,
                                        n.cores = n.cores, verbose = 1, batch = batch.factor)
  }
  
  ediff <- ediff[order(ediff[['Z']], decreasing=T),] # Reorder by Z-score
  print(sprintf("Finished %s comparison at %s", out.base, Sys.time()))  

  # Save results
  cur.fn <- file.path(home.dir, 'results', paste(MakeDateStr(), out.base, 'de.tsv', sep='.'))
  cat('#', out.base, '\n', sep=' ', file=cur.fn)
  write.table(ediff, file=cur.fn, quote=F, sep='\t', na='na', row.names=T, col.names=NA, append=T)
  
  cur.res <- list('err.mod'=err.mod, 'ediff'=ediff)
  all.de.results[[i]] <- cur.res
  names(all.de.results)[i] <- out.base
}

save.image(file.path(run.dir, "scde_de.working.RData"))
