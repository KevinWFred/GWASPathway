#!/usr/bin/env Rscript


.libPaths(c("/data/wangx53",.libPaths())) #some libs installed here

setwd("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation")

# i1 should be 1--5000 to run this procedure for 5000 simulations.
arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])

cat("Simulation with seed",i1,"\n")

library("ARTP2")

text = c("AFR","AMR","EAS","EUR","SAS")


ncase = 20000
nctrl = 20000

fam <- bim <- bed <- vector("character", 5)
studys <- vector("character", 5)
for(i in 1:5){
  p1 = paste0("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/data/refdata/",text[i],"_ref")
  file = paste0("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/result/summarydata/study_",text[i],i1,".txt.gz")
  fam[i] = paste0(p1, ".fam")
  bim[i] = paste0(p1, ".bim")
  bed[i] = paste0(p1, ".bed")
  studys[i] = file
}
pathway_file = "/data/KY_HZ_BW/pathway-T2D/SF_Kevin/pathway.txt.gz"
reference <- data.frame(fam, bim, bed, stringsAsFactors = FALSE)

# library(data.table)
# path = fread(pathway_file)
# summary  = fread(file)


## running three sARTP methods 
lambda <- 1
ncases <- list()
nctrls <- list()
ncases[[1]] <- ncase
nctrls[[1]] <- nctrl
family <- 'binomial'


# running sARTP.multiPop.SNP.union
study.list = as.list(studys)
reference.list <- split(reference, 1:nrow(reference))
lambda.list = as.list(rep(lambda,5))
ncases.list = replicate(5, ncases, FALSE)
nctrls.list = replicate(5, nctrls, FALSE)
options_union  = list()
s0 = i1+10000
for (j in 1:5){
  options_union[[j]] = list(inspect.snp.n = 2, maf = .05, HWE.p = 1e-6,
                            group.gap = 10^7, gene.R2 = .95, 
                            out.dir = "/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/result/logs",
                            nthread = 1, id.str = paste0(text[j],i1), seed = s0+j,
                            save.setup = FALSE)
}

setwd("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/result/binlogs")
res_union <- sARTP.multiPop.SNP.union_intersect(study.list, pathway_file, family,
                                                reference.list,lambda.list, 
                                                ncases.list, nctrls.list,
                                                options.list=options_union, 
                                                options.merged=list(nperm=1e4, id.str =i1))
# Removing SNPs close to marginal signals: 
# Error in update.pathway.definition(pathway, exc.snps) : 
#   All SNPs excluded in update.pathway.definition
#summaryData.setup.R line 141: ma <- find.marg.signal(sum.stat, allele.info, options)
pv.path_union = res_union$pathway.pvalue
pv.gene_union = res_union$gene.pvalue

result = list(res_union)
save(result,file=paste0("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/result/Runionres/Results_",i1,".Rdata"))
cat("SNP.union results for",i1, "are saved!\n")

# Once you have completed these 5000 simulations, 
# make sure to merge the results into a final result for the type-one errors, 
# taking into account the order of pv.gene_union during the merging process.

# tmp=data.frame(code=rep("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/Runion.R"),i1=1:5000)
# write.table(tmp,file="Runion.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#in log folder
#swarm -f /data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/Runion.swarm -g 32  --module R  --time 5:00:00 --gres=lscratch:32
