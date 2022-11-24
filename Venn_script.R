library(dplyr)
library(ggplot2)
library(reshape)
library(VennDiagram)

# load the epitranscripts genes file
epitrans_genes.df <-  read.csv("Epitrans_final.csv", header = T, sep = ",")
epitrans_genes.df <- epitrans_genes.df[,-1]

# create vector of all genes from the epitranscripts data set 
epitrans_genesv2 <- epitrans_genes.df$name
# remove the genes which are in form of integers
epitrans_genesv2 <- as.data.frame(epitrans_genesv2[3:1384])
colnames(epitrans_genesv2)[1] <- "Genes"

# load the epigenomics genes file: change the file path acc to cwd
epigen_genes.df <-  read.csv("../Epigenomics_Gull/Epigen_final.csv", header = T, sep = ",")
epigen_genes.df <- epigen_genes.df[,-1]


# create vector of all the genes from the epigenetics data set
epigen_genesv2 <- as.data.frame(epigen_genes.df$Gene.Symbol)
colnames(epigen_genesv2)[1] <- "Genes"


# get the number of different and common genes from 
#dataframes of epigenomics and epitrancriptomics study

# different genes
diff.epigen_genes<- setdiff(epigen_genesv2, epitrans_genesv2)

diff.epitrans_genes <- setdiff(epitrans_genesv2, epigen_genesv2)

# common genes
comm.epigen_genes <- intersect(epigen_genesv2, epitrans_genesv2) # 1 gene 
#comm.epitrans_genes <-intersect(epitrans_genesv2, epigen_genesv2) # 1 gene

# bind the different genes from epigenomics and epitranscriptomics
max.len <- max(c(dim(diff.epigen_genes)[1],dim(diff.epitrans_genes)[1]))

comb.diff <- data.frame(epigen.diff.genes.aftercutoff = c(diff.epigen_genes$Genes,                
                            rep(NA, max.len - dim(diff.epigen_genes)[1])),
                   epitrans.diff.genes.aftercutoff = c(diff.epitrans_genes$Genes,
                            rep(NA, max.len - dim(diff.epitrans_genes)[1])))

write.csv(comb.diff, "comb_diff_genes.csv")


max.len.all <- max(c(dim(epigen_genesv2)[1], dim(epitrans_genesv2)[1]))
comb.all <- data.frame(epigen.all.genes.aftercutoff = c(epigen_genesv2$Genes,
                        rep(NA, max.len.all -dim(epigen_genesv2)[1])),
                       epitrans.all.genes.aftercutoff = c(epitrans_genesv2$Genes,
                       rep(NA, max.len.all -dim(epitrans_genesv2)[1])))

write.csv(comb.all, "comb_all_epigen_epitrans.csv")


# venn diagram
# create pairwise Venn diagram
tiff("Venn_epigen_epitrans_genes.tiff", width = 16, height =12, res= 200, units = "in")
draw.pairwise.venn(area1=dim(epigen_genesv2)[1], area2=dim(epitrans_genesv2)[1],
                   cross.area=dim(epitrans_genesv2)[1] - dim(diff.epitrans_genes)[1],
                   category=c("Epigenomics genes","Epitranscriptomics genes"),
                   fill=c("pink", "yellow"), col="Green", lty = "dashed", scaled = T,
                   ext.text = F,
                   cat.just = list(c(1.2,-2.3), c(-2, -4)),
                   cex= 2.3,
                   cat.cex=1.5
                   )
dev.off()



