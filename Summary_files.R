#!/usr/bin/env Rscript
.libPaths(c("/data/wangx53",.libPaths())) #some libs installed here

arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
cat("Simulation with seed",i1,"\n")

library("ARTP2")
library(bigsnpr)

source("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/functions.R")
text = c("AFR","AMR","EAS","EUR","SAS")

setwd("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation")


ncase = 20000
nctrl = 20000


beds = vector("character",5)
mybeds = list()
for(i in 1:5){
  beds[i] = paste0("data/plink_select/",text[i],".bed")
  mybeds[[i]]  <- bed(beds[i])
}


# generate the summary data for warm.gene and gene.intersect
pathway_file = "/data/KY_HZ_BW/pathway-T2D/Simu/pathway.txt.gz"
pathway = read.table(gzfile(pathway_file),sep="\t",header = TRUE)

set.seed(i1)
for (j in 1:5){
  # test for the simulation data: AFR
  rds <- snp_readBed(beds[j], backingfile = tempfile())
  # Loading the data from backing files
  test <- snp_attach(rds)
  datai = test$genotypes[]
  mapi = test$map
  
  pathid = which(pathway$SNP %in% mapi$marker.ID)
  mypathway = pathway[pathid,]
  mylst <- split(1:nrow(mypathway), mypathway$Gene)
  mylst_sorted <- mylst[order(sapply(mylst, min))]
  

  control <- sampleColumnSets(datai, nctrl, mylst_sorted)
  case <- sampleColumnSets(datai, ncase, mylst_sorted)
  
  
  mydata = rbind.data.frame(control,case)
  colnames(mydata) = paste0("X", 1:ncol(mydata))
  mydata$response = c(rep(0,nctrl),rep(1,ncase))
  
  
  # fitting by logistic regression model
  res = fit_logistic_models(mydata)
  
  # matrix of summary data
  study = mapi[,c("marker.ID","chromosome","physical.pos","allele1","allele2")]
  study = cbind(study,res)
  colnames(study) <- c("SNP","Chr","Pos","EffectAllele","RefAllele","BETA","SE","P")
  
  
  file = paste0("result/summarydata0/study_",text[j],i1,".txt.gz")
  # Create a connection to a .gz file
  con <- gzfile(file, "w")
  # Save the data frame to the file
  write.table(study, file = con, row.names = FALSE, sep = "\t")
  # Close the connection
  close(con)
  cat("Summmary data for",text[j],"is saved!\n")
}

#generate swarm file
#tmp=data.frame(code=rep("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/Summary_files.R"),i1=1:5000)
#write.table(tmp,file="Summary_files.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#in log folder
#swarm -f /data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/Summary_files.swarm -g 32  --module R  --time 4:00:00 --gres=lscratch:32

