library(ggplot2)
library(dplyr)
library(reshape)
# Analysis of the peaks_ageCov.csv file to filter out the genes 
# GSE119168_Diff_peaks_ageCov.csv: downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119168


epitranscripts <-  read.csv("GSE119168_Diff_peaks_ageCov.csv", header = TRUE, sep = ",")
epitranscripts <-  as.data.frame(epitranscripts)

# get the  duplicated set of genes data frame from the original dataframe
dup.gene <-  epitranscripts[duplicated(epitranscripts$name),]

# get unique set of genes names
gene.names <- as.data.frame(unique(epitranscripts$name))
colnames(gene.names)[1] <- "genes"
# order the unique set of genes
order.gene.names <- gene.names[order(gene.names$genes),]

# sorting of the dataframe based on genes and the p-value
order.epitranscripts <- epitranscripts[order(epitranscripts$name, epitranscripts$p_value), ]

# remove the genes which are present multiple times and retain those gene among its multiplicates 
#with those with the lowest p-value
epitranscripts.2 <- data.frame()
onegeneset <- data.frame()
i <- 1
while (i < length(order.gene.names)){
  onegeneset <- filter(order.epitranscripts, order.epitranscripts$name == order.gene.names[i])
  epitranscripts.2 <- rbind(epitranscripts.2, onegeneset[1,])
  i = i +1
  if (i== length(order.gene.names)){
    onegeneset <- filter(order.epitranscripts, order.epitranscripts$name == order.gene.names[i])
    epitranscripts.2 <- rbind(epitranscripts.2, onegeneset[1,])
  }
}

write.csv(epitranscripts.2, "Epitranscripts2.csv")
# filter the genes based on p-value cutoff 
mylist= list()
gene_nr=list()
# 10e-5 is same as 10 x 10^-5 or 10 x 10-5 = 10^-4
# filtering starting from pval 10^-1
for (i in 1:10) {
  j=i+0
  print(j) # create a list of list
  mylist[[j]] <- filter(epitranscripts.2, epitranscripts.2$p_value <= 10^-j) 
  gene_nr[i] = dim(mylist[[j]])[1]
}
gene_nr= as.data.frame(gene_nr)
colnames(gene_nr) <- c("10^-1", "10^-2", "10^-3", 
                       "10^-4", "10^-5", "10^-6",
                       "10^-7", "10^-8", "10^-9", "10^-10")
rownames(gene_nr) <- c("GENE_NUMBER")


# change the dataframe shape
melt.gene_nr <- melt(gene_nr)
colnames(melt.gene_nr) <- c("pval_cutoff", "GENE_NUMBER")

# plot
tiff("Epitrans_Genes_per_pvalcutoff.tiff", width= 10, height= 8, units="in", res=400)
ggplot(melt.gene_nr, aes(x=pval_cutoff, y=GENE_NUMBER, group= 1)) +
  geom_line(color = 4) +
  geom_point(size=2) +
  theme(axis.text.x = element_text(angle = 90, size= 10)) +
  ggtitle("Number of genes after p-val cut off filtering")
dev.off()

# filter out the final data frame with p-val cutoff 10^-3
epitrans.gene_nr.final  <- mylist[[3]]
write.csv(epitrans.gene_nr.final, "Epitrans_final.csv")


# IRRELEVANT
yy <- filter(epitranscripts, epitranscripts$p_value <= 10e-4) # 1663 genes
y <-  unique(yy$name) # 1384 genes

# plots scripts adapted from:
#https://r-graph-gallery.com/about.html