#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

# Script for quickly generating a bunch of random excitatory-inhibitory donor_id pairings

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

home.dirs <- c("/home/local/users/flwu/teamawesome/", "/results/G4017_Users/team_awesome/",
               "/Users/Felix/GitHub/")
home.dir <- detectDirs(home.dirs)

data.dir <- file.path(home.dir, "datasets")
scrna.dir <- file.path(data.dir, "scRNASeq")
run.dir <- file.path(home.dir, "deepseq_awesome")
#res.dir <- file.path(home.dir, 'results')

source(file.path(run.dir, "scde_funcs.R"))

cell.info <- ReadCellInfo(file.path(scrna.dir, "columns-cells.csv"))

seed <- 182
min.cells <- 20

inhib <- c("Gad2-IRES-Cre", "Slc32a1-IRES-Cre")
excit <- c("Slc17a6-IRES-Cre")

inhib.tab <- table(cell.info[cell.info[["genotype_driver"]] %in% inhib,"donor_id"])
inhib.nms <- names(inhib.tab)[inhib.tab >= min.cells]
excit.tab <- table(cell.info[cell.info[["genotype_driver"]] %in% excit,"donor_id"])
excit.nms <- names(excit.tab)[excit.tab >= min.cells]

num.pairs <- min(length(inhib.nms), length(excit.nms))
set.seed(seed)
inhib.nms <- sample(inhib.nms, size=num.pairs)
excit.nms <- sample(excit.nms, size=num.pairs)

spec.table <- data.frame("group.col"=rep("broad_class_type", num.pairs),
                         "control"=rep("gabaergic",num.pairs),
                         "case"=rep("glutamatergic",num.pairs),
                         "filters"=paste("genotype_driver=", inhib.nms, ",", excit.nms, ";", sep=""),
                         "exclude"=rep("", num.pairs),
                         "batch.col"=rep("", num.pairs),
                         "seed"=rep("", num.pairs),
                         "out.base"=paste("glut_gaba.", inhib.nms, "_", excit.nms, sep=""))

cur.fn <- file.path(run.dir, "scde_de.rand_donor_id.specs.txt")
write.table(spec.table, file=cur.fn, quote=F, sep="\t", row.names=F)