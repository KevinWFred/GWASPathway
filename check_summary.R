#!/usr/bin/env Rscript

library(data.table)
library(ggpubr)
library("cowplot") 
qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.2,cex.axis=1.2)
  
  abline(0,1,lty=2)
  chisq <- qchisq(1-pvalue,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}

text = c("AFR","AMR","EAS","EUR","SAS")
pathway_file = "/data/KY_HZ_BW/pathway-T2D/Simu/pathway.txt.gz"
pathway = read.table(gzfile(pathway_file),sep="\t",header = TRUE)
#common snps across populations
comsnps=read.table("result/MAF05_comsnps.txt")

plotsnpZ=function(snp=fsnp)
{
  p=rep(NA,5)
  names(p)=text
  for(j in 1:length(text))
  {
    pop=text[j]
    print(pop)
    Z=rep(NA,5000)
    for (i in 1:5000)
    {
      if (i %% 500==0) cat(i,'..')
      tmp=fread(paste0("result/summarydata/study_",pop,i,".txt.gz"))
      if (i==1)
      {
        idx=which(tmp$SNP==snp)
        if (length(idx)==0) stop("snp not found")
      }
      Z[i]=tmp$BETA[idx]/tmp$SE[idx]
    }
    tmp=shapiro.test(Z)
    
    #plot1 <- recordPlot()  
    #plot2=ggqqplot(Z)
    pdf(paste0("result/checksummary/gene",gene,"_",pop,"_",snp,".pdf"))
    par(mfrow=c(1,2))
    plot(density(Z))
    #text(-1,0.2,paste0("p=",round(tmp$p.value,3)))
    qqnorm(Z)
    qqline(Z,col="red",lwd=2)
    text(-1,1.5,paste0("p=",round(tmp$p.value,3)))
    dev.off()
 
    p[j]=tmp$p.value
  }
  print(p)
  return(p)
}
gene=2551 #highest error
snps=pathway$SNP[which(pathway$Gene==gene)]
fsnp=intersect(comsnps$V1,snps)[1] #rs8133991:27109486:A:G
snps=intersect(comsnps$V1,snps) #8 common snps
fsnpp=plotsnpZ()
#snp1p=plotsnpZ(snp=snps[1])
snp2p=plotsnpZ(snp=snps[2])
snp3p=plotsnpZ(snp=snps[3])
snp4p=plotsnpZ(snp=snps[4])
snp5p=plotsnpZ(snp=snps[5])
snp6p=plotsnpZ(snp=snps[6])
snp7p=plotsnpZ(snp=snps[7])
snp8p=plotsnpZ(snp=snps[8])

gene=3275
snps=pathway$SNP[which(pathway$Gene==gene)]
fsnp=intersect(comsnps$V1,snps)[1]
snps=intersect(comsnps$V1,snps) #11 common snps
fsnpp1=plotsnpZ()
#snp1p=plotsnpZ(snp=snps[1])
snp2p1=plotsnpZ(snp=snps[2])
snp3p1=plotsnpZ(snp=snps[3])
snp4p1=plotsnpZ(snp=snps[4])
snp5p1=plotsnpZ(snp=snps[5])
snp6p1=plotsnpZ(snp=snps[6])
snp7p1=plotsnpZ(snp=snps[7])
snp8p1=plotsnpZ(snp=snps[8])

