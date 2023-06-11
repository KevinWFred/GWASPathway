#!/usr/bin/env Rscript
#create common sets of snps across 5 datsets; these snp can be used as function snps
.libPaths(c("/data/wangx53",.libPaths())) #some libs installed here

setwd("/data//KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation")

library(data.table)
check_snps=function(file1="/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/AFR.bim",
                    file2="/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/EAS.bim")
{
  dat1=as.data.frame(fread(file1))
  dat2=as.data.frame(fread(file2))
  res1=table(dat1$V2 %in% dat2$V2)
  tmp1=paste0(dat1$V1,":",dat1$V4)
  tmp2=paste0(dat2$V1,":",dat2$V4)
  res2=table(tmp1 %in% tmp2)
  #position id equal snpid, which means npid(a:b:c:d) are the same in two data
  if (any(res1!=res2)) print(paste0("snpid not consistent between ",file1," and ",file2))
}
check_snps()
check_snps(file2="/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/EUR.bim")
check_snps(file2="/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/SAS.bim")
check_snps(file2="/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/AMR.bim")
#the snpids are consistent among all dataets

check_allels=function(file1="/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/AFR.bim",
                      file2="/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/EAS.bim")
{
  dat1=as.data.frame(fread(file1))
  dat2=as.data.frame(fread(file2))
  comsnp=intersect(dat1$V2,dat2$V2)
  idx1=match(comsnp,dat1$V2)
  idx2=match(comsnp,dat2$V2)
  print(table(dat1$V5[idx1] == dat2$V5[idx2]))
  #FALSE  TRUE 
  #752  2328
  print(table(dat1$V6[idx1] == dat2$V5[idx2]))
  #FALSE  TRUE 
  #2328   752
} #the alleles are not the same

#re-generate all the data based on snpid, to make the alleles consistent
# rs184802341:9908330:A:G, the effect allele should be G
plink="/usr/local/apps/plink/1.9/plink"
regen_allels=function(prefix="AFR",outfolder="data/plink_select/")
{
  dat=as.data.frame(fread(paste0("/data/KY_HZ_BW/pathway-T2D/Simu/plink_select/",prefix,".bim")))
  #tmp=str_count(dat$V2,":")
  tmp=unlist(strsplit(dat$V2,":"))
  a1=tmp[seq(4,length(tmp),4)]
  tmp1=data.frame(snp=dat$V2,a1=a1)
  write.table(tmp1,file=paste0("result/",prefix,"_a1.txt"),row.names = F,sep=" ",quote=F)
  cmd=paste0(plink," --bfile /data/KY_HZ_BW/pathway-T2D/Simu/plink_select/",prefix," --a1-allele result/",prefix,"_a1.txt --make-bed --out ",outfolder,prefix)
  system(cmd)
  cmd=paste0(plink," --bfile ",outfolder,prefix," --keep-allele-order --recode A-transpose --out ",outfolder,prefix)
  system(cmd)
}
regen_allels()
#quick check
# allafr=as.data.frame(fread("data/plink_select/AFR.bim"))
# tmp=unlist(strsplit(allafr$V2,":"))
# tmp1=tmp[seq(4,length(tmp),4)]
# table(tmp1==allafr$V5)
# traw=as.data.frame(fread("data/plink_select/AFR.traw"))
# all(traw$COUNTED==allafr$V5) #T
regen_allels(prefix="EAS")
regen_allels(prefix="EUR")
regen_allels(prefix="SAS")  
regen_allels(prefix="AMR")

check_allels(file1="data/plink_select/AFR.bim",file2="data/plink_select/EUR.bim")

#generate data with MAF>0.05
gen_maf05=function(prefix="AFR",outfolder="data/plink_select/")
{
  cmd=paste0(plink," --bfile ",outfolder,prefix," --keep-allele-order --maf 0.05 --make-bed --out ",outfolder,prefix,"_MAF05")
  system(cmd)
}
gen_maf05()
gen_maf05(prefix="EAS")
gen_maf05(prefix="EUR")
gen_maf05(prefix="SAS")  
gen_maf05(prefix="AMR")

check_allels(file1="data/plink_select/AFR_MAF05.bim",file2="data/plink_select/EUR_MAF05.bim")

#find the common set of snps with MAF>0.05
allafr=as.data.frame(fread("data/plink_select/AFR_MAF05.bim"))
alleas=as.data.frame(fread("data/plink_select/EAS_MAF05.bim"))
alleur=as.data.frame(fread("data/plink_select/EUR_MAF05.bim"))
allsas=as.data.frame(fread("data/plink_select/SAS_MAF05.bim"))
allamr=as.data.frame(fread("data/plink_select/AMR_MAF05.bim"))
comsnps=intersect(allafr$V2,alleas$V2)
comsnps=intersect(comsnps,alleur$V2)
comsnps=intersect(comsnps,allsas$V2)
comsnps=intersect(comsnps,allamr$V2) #2441 snps
write.table(comsnps,file="result/MAF05_comsnps.txt",row.names = F,col.names = F,quote=F)

#regenerate the ref files
tmp1=fread("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/refdata/AFR_ref.bim")
tmp2=fread("data/plink_select/AFR.bim")
table(tmp1$V2==tmp2$V2)
table(tmp1$V5==tmp2$V5)
table(tmp1$V6==tmp2$V5)
regen_refallels=function(prefix="AFR",outfolder="data/refdata/")
{
  cmd=paste0(plink," --bfile /data/KY_HZ_BW/pathway-T2D/SF_Kevin/refdata/",prefix,"_ref"," --a1-allele result/",prefix,"_a1.txt --make-bed --out ",outfolder,prefix,"_ref")
  system(cmd)
}
regen_refallels(prefix="EAS")
regen_refallels(prefix="EUR")
regen_refallels(prefix="SAS")  
regen_refallels(prefix="AMR")
tmp1=fread("data/refdata/AFR_ref.bim")
tmp2=fread("data/plink_select/AFR.bim")
table(tmp1$V2==tmp2$V2)
table(tmp1$V5==tmp2$V5)
