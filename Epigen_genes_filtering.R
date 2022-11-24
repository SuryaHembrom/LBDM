# filtering of genes from the paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9327231/#MOESM1
# GSE202245
# Paper: DNA methylation and transcriptomic features are preserved throughout disease
########recurrence and chemoresistance in high grade serous ovarian cancers
library(dplyr)
library(ggplot2)
library(reshape)

supp12 <- read.csv("Supplementary Table 12_Genes near BC vs NC DMRs intersecting DEGs.csv", sep = ",", header = TRUE)
supp12 <- as.data.frame(supp12)

# genes present duplicated due to binning of the methylated regions:
supp12.dup <- supp12[duplicated(supp12$Gene.Symbol),]

# get unique set of genes names
supp12.gene.names <- as.data.frame(unique(supp12$Gene.Symbol))
colnames(supp12.gene.names)[1] <- "genes"
# order the unique set of genes
order.supp12.gene.names <- supp12.gene.names[order(supp12.gene.names$genes),]

# sorting of the dataframe based on genes and the adjusted p-value
order.supp12 <- supp12[order(supp12$Gene.Symbol, supp12$Adjusted.P.value..DESeq2.), ]

# remove the genes which are present multiple times and retain those gene among its multiplicates 
#with those with the lowest adjusted p-value
supp12.v2 <- data.frame()
supp12.onegeneset <- data.frame()
i <- 1
while (i < length(order.supp12.gene.names)){
  supp12.onegeneset <- filter(order.supp12, order.supp12$Gene.Symbol == order.supp12.gene.names[i])
  supp12.v2 <- rbind(supp12.v2, supp12.onegeneset[1,])
  i = i +1
  if(i == length(order.supp12.gene.names)){
    supp12.onegeneset <- filter(order.supp12, order.supp12$Gene.Symbol == order.supp12.gene.names[i])
    supp12.v2 <- rbind(supp12.v2, supp12.onegeneset[1,])
  }
}

write.csv(supp12.v2,"Epigen_Supp12v2.csv")

# filter the genes based on adjusted p-value cutoff 
mylist.supp12= list()
gene_nr.supp12=list()
# 10e-5 is same as 10 x 10^-5 or 10 x 10-5 = 10^-4
# filtering starting from pval  10^-1
for (i in 1:10) {
  j=i+0
  print(j) # create a list of list
  mylist.supp12[[j]] <- filter(supp12.v2, supp12.v2$Adjusted.P.value..DESeq2. <= 10^-j)
  gene_nr.supp12[i] = dim(mylist.supp12[[j]])[1]
}
gene_nr.supp12= as.data.frame(gene_nr.supp12)
colnames(gene_nr.supp12) <- c("10^-1", "10^-2", "10^-3", 
                       "10^-4", "10^-5", "10^-6",
                       "10^-7", "10^-8", "10^-9", "10^-10")
rownames(gene_nr.supp12) <- c("GENE_NUMBER")


# change the dataframe shape
melt.gene_nr.supp12 <- melt(gene_nr.supp12)
colnames(melt.gene_nr.supp12) <- c("adj_pval_cutoff", "GENE_NUMBER")

# plot
tiff("Epigenomics_Genes_per_adjpvalcutoff.tiff", width= 10, height= 8, units="in", res=400)
ggplot(melt.gene_nr.supp12, aes(x=adj_pval_cutoff, y=GENE_NUMBER, group= 1)) +
  geom_line(color = 4) +
  geom_point(size=2) +
  theme(axis.text.x = element_text(angle = 90, size= 10)) +
  ggtitle("Number of genes after adj p-val cut off filtering")
dev.off()


# filter out the final data frame with adj p-val cutoff 10^-2
#gene_nr.final  <- mylist.supp12[[2]] 

# filter out the final data frame with adj p-val cutoff 0.01
gene_nr.final <- filter(supp12.v2, supp12.v2$Adjusted.P.value..DESeq2. <= 0.01)
write.csv(gene_nr.final, "Epigen_final.csv")



# plots scripts adapted from:
#https://r-graph-gallery.com/about.html

     