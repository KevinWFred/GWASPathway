#!/usr/bin/env Rscript

.libPaths(c("/data/wangx53",.libPaths())) #some libs installed here

setwd("/data//KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation")

arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])
#causalsnpidx=as.numeric(arg[[2]]) #1-nth snp in each gene, function snp pick from common spns
causalsnpidx=1 #pick the first common snp for each gene
#beta of function snp in each pop
beta=rep(0,5)
set.seed(i1)
cat("Simulation with seed",i1,"\n")

library("ARTP2")
library(bigsnpr)

source("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/functions.R")
text = c("AFR","AMR","EAS","EUR","SAS")

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

pathway_file1 = "/data/KY_HZ_BW/pathway-T2D/SF_Kevin/pathway.txt.gz"
pathway1 = read.table(gzfile(pathway_file1),sep="\t",header = TRUE)
#common snps across populations
comsnps=read.table("result/MAF05_comsnps.txt")

for (j in 1:5){
  # test for the simulation data: AFR
  rds <- snp_readBed(beds[j], backingfile = tempfile())
  # Loading the data from backing files
  test <- snp_attach(rds)
  datai = test$genotypes[]
  mapi = test$map
  pathid = which(pathway$SNP %in% mapi$marker.ID)
  mypathway = pathway[pathid,]
  mydatai=datai[,match(mypathway$SNP,mapi$marker.ID)]
  colnames(mydatai)=mypathway$SNP
  #for each population
  mycomsnps=intersect(comsnps$V1,mypathway$SNP)
  #keep the gene order when split
  mypathway$Gene=factor(mypathway$Gene,levels=unique(mypathway$Gene))
  mylst <- split(1:nrow(mypathway), mypathway$Gene)
  mypathway$Gene=as.character(mypathway$Gene)
  #mylst_sorted <- mylst[order(sapply(mylst, min))]
  
  allcase=data.frame(matrix(NA,ncol=nrow(mypathway),nrow=ncase))
  allcontrol=data.frame(matrix(NA,ncol=nrow(mypathway),nrow=nctrl))
  colnames(allcase)=colnames(allcontrol)=mypathway$SNP
  rownames(allcase)=paste0("X",1:nrow(allcase))
  rownames(allcontrol)=paste0("X",1:nrow(allcontrol))
  allcausalsnp=NULL##keep record of all causal snps
  for (k in 1:length(mylst))
  {
    allsnp=mypathway$SNP[mylst[[k]]] #alll snps in a gene
    causalsnp=intersect(allsnp,mycomsnps)[causalsnpidx] #causal snp idx based on whole data
    if (length(causalsnp)==0) stop(paste0(names(mylst)[k]," doesn't have function snp!"))
    allcausalsnp=c(allcausalsnp,causalsnp)
    othersnp=allsnp[!allsnp %in% causalsnp] #other snps in whole data
    #allele frequencies
    f0=sum(mydatai[,causalsnp]==0,na.rm=T)/sum(!is.na(mydatai[,causalsnp]))
    f1=sum(mydatai[,causalsnp]==1,na.rm=T)/sum(!is.na(mydatai[,causalsnp]))
    f2=sum(mydatai[,causalsnp]==2,na.rm=T)/sum(!is.na(mydatai[,causalsnp]))
    ctrl_causalsnp=sample(c(0, 1, 2), size = nctrl, replace = TRUE, prob = c(f0, f1, f2))
    allcontrol[,causalsnp]=ctrl_causalsnp
    
    t1=table(ctrl_causalsnp) #counts of 0/1/2 in controls
    t11=rep(0,3)
    for (i in 1:3)
    {
      idx=which(names(t1)==i-1)
      if (length(idx)>0) t11[i]=t1[idx]
    }
    fcase=c(f0,f1*exp(beta[j]),f2*exp(2*beta[j]))
    fcase=fcase/sum(fcase)
    case_causalsnp=sample(c(0, 1, 2), size = ncase, replace = TRUE, prob = fcase)
    t2=table(case_causalsnp) #counts of 0/1/2 in cases
    t22=rep(0,3)
    for (i in 1:3)
    {
      idx=which(names(t2)==i-1)
      if (length(idx)>0) t22[i]=t2[idx]
    }
    #fill in function snp for control
    allcase[,causalsnp]=case_causalsnp
    
    #fill in other snps
    g0=which(mydatai[,causalsnp]==0) #all samples with 0
    g1=which(mydatai[,causalsnp]==1)
    g2=which(mydatai[,causalsnp]==2)
    t3=t11+t22 #counts of 0/1/2 in controls + cases #assume t3 has 0/1/2
    #genotype 0
    if (length(g0)<t3[1]) stop(paste0("not enough 0 for ",causalsnp))
    othersample0idx=sample(g0,t3[1]) #sampling 0 gt in case and control
    ctrl_othersample0idx=othersample0idx[1:t11[1]] #for control
    case_othersample0idx=othersample0idx[!othersample0idx %in% ctrl_othersample0idx]
    #genotype 1
    if (length(g1)<t3[2]) stop(paste0("not enough 1 for ",causalsnp))
    othersample1idx=sample(g1,t3[2])
    ctrl_othersample1idx=othersample1idx[1:t11[2]]
    case_othersample1idx=othersample1idx[!othersample1idx %in% ctrl_othersample1idx]
    #genotype 2
    if (length(g2)<t3[3]) stop(paste0("not enough 2 for ",causalsnp))
    othersample2idx=sample(g2,t3[3])
    ctrl_othersample2idx=othersample2idx[1:t11[3]]
    case_othersample2idx=othersample2idx[!othersample2idx %in% ctrl_othersample2idx]
    #fill in the data
    allcontrol[,othersnp]=mydatai[c(ctrl_othersample0idx,ctrl_othersample1idx,ctrl_othersample2idx),othersnp]
    allcase[,othersnp]=mydatai[c(case_othersample0idx,case_othersample1idx,case_othersample2idx),othersnp]
  }
  print(paste0(Sys.time()," genotype matrix is ready"))
  #row: sample, column:snps
  control <- allcontrol
  case <- allcase
  
  mydata = rbind.data.frame(control,case)
  colnames(mydata) = paste0("X", 1:ncol(mydata))
  mydata$response = c(rep(0,nctrl),rep(1,ncase))
  
  # fitting by logistic regression model
  res = fit_logistic_models(mydata)
  
  # matrix of summary data
  study = mapi[,c("marker.ID","chromosome","physical.pos","allele1","allele2")]
  idx=match(colnames(control),mypathway$SNP)
  study = cbind(study,res[idx,])
  colnames(study) <- c("SNP","Chr","Pos","EffectAllele","RefAllele","BETA","SE","P")
  
  file = paste0("result/summarydata/study_",text[j],i1,".txt.gz")
  # Create a connection to a .gz file
  con <- gzfile(file, "w")
  # Save the data frame to the file
  write.table(study, file = con, row.names = FALSE, sep = "\t")
  # Close the connection
  close(con)
  print(Sys.time())
  cat("Summmary data for",text[j],"is saved!\n")
}
#generate swarm file
#tmp=data.frame(code=rep("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/generate_summary.R"),i1=1:5000)
#write.table(tmp,file="generate_summary.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#in log folder
#swarm -f /data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/generate_summary.swarm -g 32  --module R  --time 4:00:00 --gres=lscratch:32
#check results
# tmp=list.files("result/summarydata/","study_AFR\\w+.txt.gz")
# length(tmp) #4895
# tmp=list.files("result/summarydata/","study_EAS\\w+.txt.gz")
# length(tmp)
# tmp=list.files("result/summarydata/","study_EUR\\w+.txt.gz")
# length(tmp)
# tmp=list.files("result/summarydata/","study_SAS\\w+.txt.gz")
# length(tmp)
# tmp1=gsub("study_SAS","",tmp)
# tmp1=as.numeric(gsub(".txt.gz","",tmp1,fixed=T))
# tmp2=1:5000
# tmp3=tmp2[!tmp2 %in% tmp1]
# tmp4=data.frame(code=rep("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/generate_summary.R"),i1=tmp3)
# write.table(tmp4,file="generate_summary_missing.swarm",row.names = F,col.names = F,sep=" ",quote=F)
#swarm -f /data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/generate_summary_missing.swarm -g 64  --module R  --time 5:00:00 --gres=lscratch:32
