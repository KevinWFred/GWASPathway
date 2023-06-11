#!/usr/bin/env Rscript

# source("/data/KY_HZ_BW/pathway-T2D/Simu/functions.R")

setwd("/data/KY_HZ_BW/pathway-T2D/SF_Kevin/power_simulation/result/Runionres/")

# Select all files with the extension ".txt" in the current directory
files <- list.files(pattern = "*.Rdata", full.names = TRUE)
# files <- list.files(path = ".", pattern = "*.Rdata", full.names = TRUE)
# Extract the key number from the file names
files1 <- files[grep("Results", basename(files))]
key_numbers <- as.numeric(gsub("\\D", "", basename(files1)))

# Sort the files by the key numbers
files1 <- files1[order(key_numbers)]
seeds = key_numbers[order(key_numbers)]

n.files = length(files1)
seedvec <- Pv.path <- rep(0,n.files)
  
load(files1[1])
genes = result[[1]]$gene.pvalue$Gene
genes = genes[order(genes)]
ngene = length(genes)
Pv.gene = matrix(NA, nrow = ngene, ncol = n.files)
rownames(Pv.gene) = genes



for(i in 1:n.files){
  load(files1[i])
  if (i %%100==0){
    cat("File",i,"is loaded!!\n")
  }
  
  Pv.path[i] = result[[1]]$pathway.pvalue
  df = result[[1]]$gene.pvalue
  id = match(df$Gene,genes)
  Pv.gene[id,i] = df$Pvalue
  
  seedvec[i] = seeds[i]
  rm(result)
}

# print(length(seedvec)==2000)
res.path = mean(Pv.path<=5E-02, na.rm = TRUE) #0.081
print(res.path)
# quantile(Pv.path,c(0,0.05,0.1,0.25,0.5,1))
# 0%          5%         10%         25%         50%        100% 
#   0.000149985 0.025594941 0.067438256 0.242300770 0.586066393 1.000000000  

res.gene = rep(0,ngene)
names(res.gene) = genes
res.gene = rowMeans(Pv.gene <=5E-02, na.rm = TRUE)
quantile(res.gene,c(0,0.25,0.5,0.75,0.9,0.95,1))
# 0%     25%     50%     75%     90%     95%    100% 
# 0.00000 0.00195 0.00400 0.01375 0.10462 0.22809 0.52580

which.max(res.gene) #29
res.gene[29] #2551

which.min(res.gene) #37
res.gene[37] #3275
save(Pv.path,Pv.gene,file = "../NULL_Union.Rdata")

res = c(res.path,res.gene)
names(res)[1] = "Pathway"
df_res = data.frame(union = res)

write.csv(df_res, file = "../NULL_5k_union.csv")




