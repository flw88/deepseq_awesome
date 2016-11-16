#!/usr/bin/env Rscript

# Convert FPKM to counts using some sort of heuristic or back-calculation

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

library(monocle)

home.dirs <- c("/home/local/users/flwu/teamawesome/", "/results/G4017_Users/team_awesome/",
               "/Users/Felix/GitHub/")
home.dir <- detectDirs(home.dirs)

scrna.dir <- file.path(home.dir, "datasets", "scRNASeq")

fpkm.df <- read.csv(file.path(scrna.dir, "fpkm_table.csv"), header=T, row.names=1)
names(fpkm.df) <- gsub("^X", "", x=names(fpkm.df)) # Remove preceding 'X' that R inserts in col names
gene.info <- read.csv(file.path(scrna.dir, "rows-genes.csv"), header=T, row.names=1, quote='"')
cell.info <- read.csv(file.path(scrna.dir, "columns-cells.csv"), header=T, row.names=1)

# Check that all the identifiers correspond exactly
stopifnot(all(names(fpkm.df) == rownames(cell.info)))
stopifnot(all(rownames(fpkm.df) == rownames(gene.info)))

# pd <- new("AnnotatedDataFrame", data=cell.info)
# fd <- new("AnnotatedDataFrame", data=gene.info)
# 
# cds <- newCellDataSet(as.matrix(fpkm.df), phenoData=pd, featureData=fd)
# 
# # Convert to absolute values using 'monocle' pkg
# use.fn <- file.path(scrna.dir, "count_table.monocle.csv")
# if(!file.exists(use.fn)){
#   count.matrix <- relative2abs(cds)
#   
#   # Write new count matrix to file
#   new.col.nm <- c('"gene_id \ rnaseq_profile_id"', colnames(count.matrix))
#   count.matrix <- cbind(rownames(count.matrix), count.matrix)
#   colnames(count.matrix) <- new.col.nm
#   write.table(count.matrix, file=use.fn, sep=",", quote = F, row.names = F, col.names = T)
# } else {
#   print(sprintf("%s exists, skipping..."))
# }

# Convert to absolute values simply by taking the ceiling
use.fn <- file.path(scrna.dir, "count_table.ceiling.csv")
if(!file.exists(use.fn)){
  count.matrix <- ceiling(as.matrix(fpkm.df))
  
  # Write new count matrix to file
  new.col.nm <- c('"gene_id \ rnaseq_profile_id"', colnames(count.matrix))
  count.matrix <- cbind(rownames(count.matrix), count.matrix)
  colnames(count.matrix) <- new.col.nm
  print(sprintf("Writing %s ..."))
  write.table(count.matrix, file=use.fn, sep=",", quote = F, row.names = F, col.names = T)
  print("Finished.")
} else {
  print(sprintf("%s exists, skipping."))
}

