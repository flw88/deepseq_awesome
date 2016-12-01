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

library(methods)
library(scde)
library(digest)

print(sprintf("Started scde_de.R at %s", Sys.time()))

home.dirs <- c("/home/local/users/flwu/teamawesome/", "/results/G4017_Users/team_awesome/",
               "/Users/Felix/GitHub/")
home.dir <- detectDirs(home.dirs)

data.dir <- file.path(home.dir, "datasets")
scrna.dir <- file.path(data.dir, "scRNASeq")
run.dir <- file.path(home.dir, "deepseq_awesome")
res.dir <- file.path(home.dir, 'results')

source(file.path(run.dir, "scde_funcs.R"))

count.fn <- "count_table.fpkm2count.csv"
count.df <- read.csv(file.path(scrna.dir, count.fn), header=T, row.names=1)
gene.info <- read.csv(file.path(scrna.dir, "rows-genes.csv"), row.names=1, header=T, quote='"')
cell.info <- ReadCellInfo(file.path(scrna.dir, "columns-cells.csv"))

spec.table <- read.delim(file.path(run.dir, "scde.Snap25.batch.specs.txt"), header=T, sep='\t', comment.char='#', as.is=T, colClasses='character')

n.cores <- 16

#Change id names
names(count.df) <- gsub("^X", replacement = '', names(count.df))

#Order info dfs
stopifnot(all(rownames(count.df) %in% rownames(gene.info)))
gene.info <- gene.info[rownames(count.df),]
stopifnot(all(names(count.df) %in% rownames(cell.info)))
cell.info <- cell.info[names(count.df),]

for(i in 1:nrow(spec.table)){
  #Read in, parse spec table info
  spec.row <- spec.table[i,]
  out.base <- spec.row[['out.base']]
  
  #Memoizing to speed things so we don't repeat comparisons needlessly
  dig <- digest(c(count.fn, spec.row))
  res.fn <- paste(dig, 'RData',sep='.')
  res.path <- file.path(res.dir, res.fn)
  if(file.exists(res.path)){
    print(sprintf("%s exists. No need to do comparison %s", res.fn, out.base))
    print(sprintf("Results can be loaded from:"))
    print(sprintf("%s\n", res.path))
  }else{
    # Run the comparison
    group.col <- spec.row[['group.col']]
    class.pair <- c(spec.row[['case']], spec.row[['control']])
    names(class.pair) <- c('case', 'control')
    
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
    
    if((spec.row[['seed']] != "") && (spec.row[['seed']] != ".") && !is.null(spec.row[['seed']])){
      print("Seed specified. Randomly selecting a cells to balance our comparisons.")
      use.lgcl <- list('case'=cur.cell.info[[group.col]] %in% class.pair['case'],
                       'control'=cur.cell.info[[group.col]] %in% class.pair['control'])
      use.n <- lapply(use.lgcl, sum)
      if(use.n[['case']] >= use.n[['control']]){
        reduce.nm <- 'case'
        keep.nm <- 'control'
      } else if (use.n[['case']] < use.n[['control']]){
        reduce.nm <- 'control'
        keep.nm <- 'case'
      }
      set.seed(as.integer(spec.row[['seed']]))
      reduce.cells <- sample(which(use.lgcl[[reduce.nm]]), size = use.n[[keep.nm]], replace=FALSE)
      keep.cells <- sort(c(reduce.cells, which(use.lgcl[[keep.nm]])))
      
      cur.cell.info <- cur.cell.info[keep.cells,]
      cur.count.m <- cur.count.m[,keep.cells]
    }
    
    batch.col <- spec.row[['batch.col']]  # which variables do we want to use as batch vars?
    if(batch.col == '' || batch.col == '.'){
      batch.factor <- NULL
    } else {
      batch.factor <- factor(cur.cell.info[[batch.col]], levels=unique(cur.cell.info[[batch.col]]))
      names(batch.factor) <- rownames(cur.cell.info)
    }
    
    print(sprintf("Started %s comparison at %s", out.base, Sys.time()))
    
    # Get factor of DE group
    group.factor <- factor(cur.cell.info[[group.col]], levels=class.pair)
    names(group.factor) <- rownames(cur.cell.info)
    
    # Clean count frame
    cur.clean.m <- clean.counts(cur.count.m, min.lib.size=1000, min.reads = 1, min.detected = 1)
    
    # Calculate error model
    err.mod <- scde.error.models(counts = cur.clean.m, groups = group.factor, n.cores = n.cores,
                                 threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
    
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
    
    ediff[['p_value']] <- 2*pnorm(abs(ediff[['Z']]),lower.tail=F) # 2-tailed p-value
    ediff[['FDR']] <- 2*pnorm(abs(ediff[['cZ']]),lower.tail=F) # FDR-corrected
    ediff[['gene_symbol']] <- gene.info[row.names(ediff), 'gene_symbol']
    ediff[['ranking.stat']] <- computeRankingStat(ediff)
    
    ediff <- ediff[order(ediff[['p_value']], decreasing=T),] # Reorder by p.value
    print(sprintf("Finished %s comparison at %s", out.base, Sys.time()))
    
    # Save results
    cur.fn <- file.path(home.dir, 'results', paste(MakeDateStr(), out.base, 'de.tsv', sep='.'))
    cat('#', out.base, '\n', sep=' ', file=cur.fn)
    write.table(ediff, file=cur.fn, quote=F, sep='\t', na='na', row.names=T, col.names=NA, append=T)
    
    cur.res <- list('err.mod'=err.mod, 'ediff'=ediff, 'cell_ids'=cur.cell.info[['rnaseq_profile_id']], 'spec.row'=spec.row)
    
    print(sprintf("Saving %s comparison results to %s", out.base, res.path))
    save(cur.res, file=res.path)
  }
}
