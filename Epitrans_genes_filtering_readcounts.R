library(dplyr)
library(DESeq2)
library(edgeR)
library(reshape)
library(plyr)
library(tidyverse)
library(arsenal)
library(statmod)
library(plotly)
library(ggfortify)
library(caret)
library(tidymodels)
library(e1071)
library(genefilter)
library(Boruta)
library(randomForest)
library(ggtext)



# the dataset is from PE reads submitted in NCBI GEO dataset as GSE119168
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119168

# download the homo sapiens annotation gtf file from Ensembl: https://www.ensembl.org/Homo_sapiens/Info/Index
# 
epitrans.counts <-  gzfile("GSE119168_readCounts.tsv.gz", "rt")
epitrans.counts.df <-  read.table(epitrans.counts, header = TRUE, sep= "\t")

# this data frame contains the number of reads mapping to exomic regions of gene binned at 50 bp 
summary(epitrans.counts.df)

# remove the binning and combine the read counts for a given gene
#1. modify the rownames of the dataframe with bins size in the rownames
    ## create unique list of gene names
x <-  list(c(gsub(",.*", "", rownames(epitrans.counts.df))))
names(x) <- "genes"
epitrans.counts.df <- cbind(epitrans.counts.df, x)
# sum the read counts per sample based on the genes
epitrans.counts.bind.df <- ddply(epitrans.counts.df, "genes", numcolwise(sum))

# to get the genes used in this sequencing 
x.df <- as.data.frame(epitrans.counts.bind.df$genes)
colnames(x.df) <- "GENES" 

write.table(x.df$GENES, "genelist_dataset.txt")

# rename the rownames based on genes from genes column 
for (i in 1:length(epitrans.counts.bind.df$genes)){
    rownames(epitrans.counts.bind.df)[i] <- epitrans.counts.bind.df$genes[i]
  }

epitrans.counts.bind.df <- epitrans.counts.bind.df[,-1]
# to get the per million reads of each sample for all the genes 
epitrans.sample.counts <- as.data.frame(colSums(epitrans.counts.bind.df)/1e6)
colnames(epitrans.sample.counts)[1] <- c("Total_Reads_in_millions")

# load the gff file 

gff <- read.table("Homo_sapiens.GRCh38.108.chr.gff3.V2.txt")
colnames(gff) <- c("chr_nr", "ensembl", "feature", "start", "end", "strand", "gene_descr", "gene_name")

# filter the gff dataframe based on genes found in the dataset

name.filtered.df  <-  c("feature", "start", "end", "hgnc_genename")
filtered.gff <- data.frame(matrix(nrow = 0, ncol = length(name.filtered.df)))
colnames(filtered.gff) <- name.filtered.df

filtered.gff <- filter(gff, gene_name %in% x.df$GENES)

filter.epitrans.counts <- with(epitrans.counts.bind.df, epitrans.counts.bind.df[rownames(epitrans.counts.bind.df)
                                                            %in% filtered.gff$gene_name,])
# gene length 
all.gene.length.gff <- as.data.frame(filtered.gff$end- filtered.gff$start)
names(all.gene.length.gff) <-  "gene_length"

filtered.gff <- cbind(filtered.gff, all.gene.length.gff)
colnames(filtered.gff)[9] <- "gene_length"

# create DESeq object 
condition.disease <- read.csv("sorted_condition_samplev2.txt", header = FALSE, row.names = 1, sep= ",")
names(condition.disease) <- "condition"

design.disease <- factor(c(condition.disease$condition))

dds <- DESeqDataSetFromMatrix(countData = filter.epitrans.counts,
                              tidy= FALSE,
                              ignoreRank = FALSE,
                              DataFrame(design.disease),
                               ~design.disease)
View(counts(dds))
# Normalisation factors for each sample
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# normalised counts matrix 
normalised_counts.deseq <- as.data.frame(counts(dds, normalized = TRUE))
pdf("normalisedcounts.pdf", width= 8.27, height = 6.5)
boxplot(normalised_counts.deseq, las =2, cex.axis=0.6, ylab = "CPM")
dev.off()

# log 
log.counts.deseq <- log2(normalised_counts.deseq+1)
pdf("normalised_logcounts.pdf")
boxplot(log.counts.deseq, las =2, cex.axis=0.6, ylab = "CPM (log2)")
dev.off()
#### get the gene length for set of genes for read counts
gene.filter.epitrans.counts <- data.frame(rownames(filter.epitrans.counts))
names(gene.filter.epitrans.counts) <- "genes"

write.table(gene.filter.epitrans.counts, "gene_filter_epitrans_counts.txt")

# order the filtered.gff file
filtered2.gff <- filtered.gff[order(filtered.gff$gene_name), ]

filtered3.gff <- filtered2.gff[, c('gene_name', 'gene_length')]
write.csv(filtered3.gff, "filtered3.gff")

# shell scripting:  to generate file with unique values with highest gene length for a given gene
# while read line; 
#do c=`grep -Ei -w "$line" filtered4.gff | sort -k2 -r |awk 'NR== 1'`; 
#echo $c >> filtered_final.txt; done < filtered3.genelist.txt

names.final <- c("gene_name", "gene_length")

filtered_final <- read.table("filtered_final.txt")
names(filtered_final) <- names.final

final.genelist <- data.frame(matrix(nrow=0, ncol=2))
names(final.genelist) <- names.final

# filter out the genes which are present multiple times in the .gff file due to their various forms present 
# as lnc-gene and pseudogene, so taking the gene feature which is more reliable 
for (i in 1:dim(gene.filter.epitrans.counts)[1]){
  final.genelist[i,] <- filtered3.gff[filtered3.gff$gene_name %in% gene.filter.epitrans.counts$genes[i], ]
}

# compare the order of genes in filter.epitrans.counts and final.genelist from the gff file
genenames <- data.frame(final.genelist$gene_name)
comparedf(gene.filter.epitrans.counts, genenames )
all.equal(gene.filter.epitrans.counts, genenames)
comparedf(filtered_final, final.genelist) # both the files produced from the R and shell are same


# get the TPM (transcripts per million) instead of RPKM or FPKM as RPKM is only for single reads and 
# TPM is more relevant than FPKM as sum of all TPM for all samples is same 
# see https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/ for more details)
# the dataset here is paired end sequenced data so TPM is calculated

#1. calculate the RPK (reads per kilobase) for each gene 
rpk <- (filter.epitrans.counts*(10^3)/(final.genelist$gene_length))
#2. get the per million scaling factor
pm.scalingfactor <- as.data.frame(colSums(rpk)/1e6)
names(pm.scalingfactor) <- "per_million_scaling_factor"
#3. calculate the TPM 
TPM <- (rpk)/(pm.scalingfactor$per_million_scaling_factor)
# check if the TPM of all the samples are the same
TPM_persample <- as.data.frame(colSums(TPM))

# getTMM 
# create a DGEList based on the grouping of disease condition
rpk.norm <- DGEList(counts = rpk, group =  design.disease)
# calcNorm : normalize the library size based on scaling factors that minimise the 
# log2 fold changes between samples for genes. 
# see edgeR manual for more details https://bioconductor.org/packages/release/bioc/html/edgeR.html

# Trimmed mean of M values (TMM) calculation for each pair of samples:
## is done for normalisation for removing composition biases between libraries:
rpk.norm <- calcNormFactors(rpk.norm)

# convert the normalised rpk into counts per million 
counts_GeTMM <- as.data.frame(cpm(rpk.norm))
pdf("counts_GeTMM.pdf", height= 7, width = 9)
boxplot(counts_GeTMM, ylim=c(0,4000), las =2, cex.axis=0.6, ylab="counts TMM")
dev.off()
# Log geTMM
log.counts_GeTMM <- log2(counts_GeTMM+1)
pdf("counts_GeTMMlog.pdf")
boxplot(log.counts_GeTMM, las =2, cex.axis=0.6, ylab="counts TMM(log2)")
dev.off()
 ###################################################

# Differential Gene Expression Analysis 
# dispersion values and gene wise negative binomial GLM with quasi likelihood tests
design <- model.matrix(~ design.disease)

data.dispersion <- estimateDisp(rpk.norm, design, robust = TRUE)
pdf("BCVplot.pdf")
plotBCV(data.dispersion)
dev.off()

# gene expression analysis from edgeR package
fit <- glmQLFit(data.dispersion, design, robust=TRUE)
pdf("QLDisp.pdf")
plotQLDisp(fit)
dev.off()

# the differentially expressed genes test
de.genes <-  glmQLFTest(fit)
# top down-regulated and up-regulated expressed genes in psoriasis sorted by log-fold change and p value cutoff 0.001
summary(decideTests(de.genes))
topTag <- topTags(de.genes,  sort.by = "logFC", p.value = 0.001) # down and up regulated genes 
de.genes.df <- as.data.frame(decideTests(de.genes))

# insignificant genes
insignificant.genes <- rownames(de.genes.df)[which(de.genes.df$design.diseaseomental_tumor ==0)]

# down-regulated genes list 
downregulated.genes <- rownames(de.genes.df)[which(de.genes.df$design.diseaseomental_tumor== -1)]

# upregulated genes list 
upregulated.genes <- rownames(de.genes.df)[which(de.genes.df$design.diseaseomental_tumor== 1)]

# filtering of 
# for filtering out the non significant genes 
insig.geneslist <- list(insignificant.genes)
insig.geneslist <- unlist(insig.geneslist)

# transpose the normalised TMM counts
transpose.counts_GeTMM <- as.data.frame(t(counts_GeTMM))
# filter out the insignificant genes from the normalized TMM transposed counts
transpose.counts_GeTMM.filt <- transpose.counts_GeTMM[, !(names(transpose.counts_GeTMM) %in% insig.geneslist)]

# pca for upregulated and downregulated genes
pca <- prcomp(transpose.counts_GeTMM.filt, scale. = TRUE)
summary(pca)
screeplot(pca)


# bind the data frames for PCA analysis
transpose.counts_GeTMM.filt.discond <- cbind(transpose.counts_GeTMM.filt, condition.disease)

#  plot the PCA of up-regulated and down-regulated genes
plot.pca <- autoplot(pca, data = transpose.counts_GeTMM.filt.discond, colour = 'condition') + 
  labs(title = "Upregulated and Downregulated genes")
ggplotly(plot.pca)


#######################################  
# get the T test for the genes which are significant between pairs of samples comparison

genes.ttest.regmodels <- rowttests(as.matrix(t(transpose.counts_GeTMM.filt)), design.disease)
ttest.raw <- genes.ttest.regmodels$p.value

# multiple hypothesis testing correction and adjustment of individual p values using Benjamini-hochberg and bonferroni methods
ttest.bh <-p.adjust(genes.ttest.regmodels$p.value, "BH", n= length(genes.ttest.regmodels$p.value))
ttest.bf <- p.adjust(genes.ttest.regmodels$p.value, "bonferroni", n = length(genes.ttest.regmodels$p.value))
res.pvals <- as.data.frame(cbind(ttest.raw, ttest.bf, ttest.bh))

# it will be good to consider the BH method and not the bonferroni that works based on family wise error rate 
# therefore it is better to use the BH p-val correction that adjusts the p-val

# ttest raw: number of genes selected 
sel.genes.raw <- which(res.pvals$ttest.raw<0.05)# 2638 genes retained
# ttest bf: number of genes selected 
sel.genes.bf <- which(res.pvals$ttest.bf<0.05) # 176 genes retained
# ttest bh: number of genes selected 
sel.genes.bh <- which(res.pvals$ttest.bh < 0.05) # 2617 genes retained
# retain the dataset with benjamini-hochberg p-val adjustment
ttbh.transpose.counts_GeTMM.filt <- transpose.counts_GeTMM.filt[ ,sel.genes.bh]


# Filter the genes to retain only the important ones and the downregulated genes for further analysis 
# filtering out of genes whose expression not > 0.0 expression in atleast 20% of the samples
ff1 <- filterfun(pOverA(p=0.20,A=0.0))
# scaling of data is not applied to the random forest
gf1 <- genefilter(ttbh.transpose.counts_GeTMM.filt,ff1) 
reduced.gset <- ttbh.transpose.counts_GeTMM.filt[gf1,]
#reduced.gset <- log2(ttbh.transpose.counts_GeTMM.filt[gf1,])

reduced.gset.dc <- cbind(reduced.gset, condition.disease)
# rename the colnames of genes as Random Forest cannot find it  
names(reduced.gset.dc) <- gsub("\\-", "\\_", names(reduced.gset.dc))

# Boruta algorithm for feature selection:
#adapted from: https://finnstats.com/index.php/2021/05/03/random-forest/
set.seed(11)
boruta <- Boruta( as.factor(condition) ~ ., data = reduced.gset.dc, doTrace = 2, maxRuns = 3000)
print(boruta)
tent.boruta <- TentativeRoughFix(boruta)
print(tent.boruta)

tent.boruta.attr <- getSelectedAttributes(tent.boruta, withTentative = F)

# plot the boruta 
pdf("boruta.pdf", width= 35, height = 20)
plot(boruta, las = 2, cex.axis = 0.5, xlab = NULL)
dev.off()
# plot importance of boruta
plotImpHistory(boruta)


# Random forest : supervised learning
###Data partition into training and test 
set.seed(20)
rf.sample <- initial_split(reduced.gset.dc, prop=0.8) # 80-20% 
rf.train <- training(rf.sample)
rf.test <- testing(rf.sample)

# Boruta on training data 
train.boruta <- Boruta(as.factor(condition) ~ ., data = rf.train, doTrace = 2, maxRuns = 3000)
print(train.boruta)

tent.train.boruta <- TentativeRoughFix(train.boruta)
print(tent.train.boruta)

train.boruta.attr<- getSelectedAttributes(train.boruta, withTentative = F)
df.train.boruta.attr <- as.data.frame(train.boruta.attr)

## THIS METHOD INCREASES THE NUMBER OF RUNS NEEDED TO SELECT THE IMPORTANT FEATURES AND HENCE MAKES THE 
## COMPUTATION TIME COMPLEX AND PRODUCES LESS IMPORTANT FEATURES AS THE SIZE OF THE DATASET IS LOW FOR 
## CREATING THE SHADOW FEATURES. 
## DECREASE IN DATASET SIZE RETURNS LESSER NUMBER OF IMPORTANT FEATURES
# WHEREAS THE ACCURACY WILL REMAIN INCREASED WITH LESSER NUMBER OF IMPORTANT FEATURES

# get features for random Forest from Boruta run on train data
trf.train.boruta_genes <- getNonRejectedFormula(tent.train.boruta)

rf.trf.train.boruta <- randomForest(trf.train.boruta_genes,data = rf.train, 
                                    importance= TRUE, mtry= length(train.boruta.attr), ntree= 10000) 
rf.trf.train.boruta



# OOB predictor try at each split for random forest
n.pred.rf <- ncol(reduced.gset.dc) -1 # removing the condition column and with -1 less predictor

# random forest on training data setwithout boruta algorithm feature selection
rf1 <- randomForest(as.factor(condition) ~., data = rf.train, 
                    importance= TRUE, mtry= n.pred.rf, ntree= 10000)
rf1 # OOB error rate is 5%
# predict the rf model in the test data
yhat.rf1 <- predict(rf1, rf.test, type = 'response')
# confusion matrix
confusionMatrix(yhat.rf1, as.factor(rf.test$condition))

# importance of features selected  
imp.rf1 <- as.data.frame(importance(rf1, type=1, scale = F)) # get the mean decrease in accuracy
# sort the dataframe according to the decrease in the mean decrease accuracy
impsort.rf1 <-  imp.rf1 %>% arrange(desc(MeanDecreaseAccuracy))

# plot the importance of the features selected 
varImpPlot(rf1, sort = TRUE, cex= 0.75)

# get features for random Forest from Boruta run
boruta_genes <- getNonRejectedFormula(tent.boruta)

# random forest on training dataset with boruta algoritm run on full dataset
rf2 <- randomForest(boruta_genes,data = rf.train, importance= TRUE, mtry= length(tent.boruta.attr), ntree= 10000) 
# 280 as the genes considered important by the boruta algorithm
rf2 # OOB error rate is 0% 
yhat.rf2 <- predict(rf2, rf.test, type = 'response')
confusionMatrix(yhat.rf2, as.factor(rf.test$condition))

# importance of features selected 
imp.rf2 <- as.data.frame(importance(rf2, type=1, scale = F)) # get the mean decrease in accuracy
# sort the dataframe  according to the mean decrease accuracy
impsort.rf2 <-  imp.rf2 %>% arrange(desc(MeanDecreaseAccuracy))

# plot the importance 
pdf("rf2_varimp.pdf")
varImpPlot(rf2, sort = TRUE, cex= 0.75)
dev.off()

# filter the genes which are important based on mean decrease accuracy
filtg.rf2 <- subset(impsort.rf2, impsort.rf2$MeanDecreaseAccuracy > 0.0) 


# RENAME THE GENES BACK AGAIN 
rownames(filtg.rf2) <- gsub("\\_", "\\-", rownames(filtg.rf2))

# FINAL GENE LIST AFTER BORUTA+ TENTATIVE FIX +RANDOM FOREST 
filtg.rf2.list <- as.data.frame(rownames(filtg.rf2))
names(filtg.rf2.list) <- c("genes")
write.csv(filtg.rf2.list, "FINAL_GENES_EPITRANS.csv")


# full dataset randomforest 
rf.full <- randomForest(as.factor(condition) ~., data = reduced.gset.dc,
                        importance= TRUE, mtry= n.pred.rf,ntree= 10000)
rf.full 

varImpPlot(rf.full, sort= TRUE, cex = 0.75)
# abort full dataset randomforest as the partition of data is important to know the
# prediction levels of the model for a new dataset


# HEATMAP 
# rename the genes of the reduced.gset.dc 
names(reduced.gset.dc) <- gsub("\\_", "\\-", names(reduced.gset.dc))

# filter out only the genes selected from random forest rf2 to obtain a dataframe
rf2.final.df <- reduced.gset.dc[, (names(reduced.gset.dc) %in% rownames(filtg.rf2))]
# add the disease condition as a new column 
rf2.final.df.dc <- cbind(rf2.final.df, condition.disease)


# COMPLETE DATSET HEATMAP #

# arrange the dataframe based on normal fallopian tube samples followed by omental tumor samples
samples.arrange <- c('OvCa1850.input', 'OvCa1850.IP', 'OvCa1917.input',
         'OvCa1917.IP', 'OvCa2005.input', 'OvCa2005.IP',
         'OvCa2013.input', 'OvCa2013.IP', 
         'OvCa2053.input', 'OvCa2053.IP', 'OvCa2064.input', 
         'OvCa2064.IP', 'OvCa2072.input', 'OvCa2072.IP', 
         'OvCa2380.IP','OvCa2380.input', 'OvCa2343.IP',
           'OvCa2343.input', 'OvCa2270.IP', 'OvCa2270.input',
           'OvCa2261.IP', 'OvCa2261.input', 'OvCa2221.IP', 
           'OvCa2221.input','OvCa2186.IP', 'OvCa2186.input')

arr.rf2.final.df.dc <- rf2.final.df.dc %>% arrange(factor(rownames(rf2.final.df.dc), 
                   levels = samples.arrange))

# ## melt the dataframe
arr.rf2.final.df.dc  <- arr.rf2.final.df.dc %>% select(-condition)
melt.arr.rf2.final.df.dc <- melt(arr.rf2.final.df.dc)
melt.arr.rf2.final.df.dc$samples <- rep(row.names(arr.rf2.final.df.dc),length(arr.rf2.final.df.dc))

# normal fallopian tube samples 
nft = c('OvCa1850.input', 'OvCa1850.IP', 'OvCa1917.input',
  'OvCa1917.IP', 'OvCa2005.input', 'OvCa2005.IP',
  'OvCa2013.input', 'OvCa2013.IP', 
  'OvCa2053.input', 'OvCa2053.IP', 'OvCa2064.input', 
  'OvCa2064.IP', 'OvCa2072.input', 'OvCa2072.IP')

myfilter <- ifelse(melt.arr.rf2.final.df.dc$samples %in% nft, "blue", "red")
# blue - normal fallopian tube
# red - omental tumor

pdf("heatmap.pdf", width = 12, height = 35)
ggplot(melt.arr.rf2.final.df.dc, aes(samples, variable)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low="white", high="black")  +
  theme(axis.text.x = element_markdown(angle = 90, size = 8, color = myfilter), 
        axis.text.y = element_text(size=9)) +
  labs(y= "Genes", x="Samples") + theme(legend.position = "none") 
dev.off()

#PLOT HEATMAP FOR GENES COMMON TO TRAINING SET BORUTA+ TENTATIVE FIX + RANDOM FOREST and FULL DATA SET + 
# TENTATIVE FIX + RANDOM FOREST 
# intersect genes from both boruta runs
genes.intersect <- intersect(train.boruta.attr, tent.boruta.attr) 
# filter out only the genes selected from random forest rf2 to obtain a dataframe
intersect.rfs.final.df <- reduced.gset.dc[, (names(reduced.gset.dc) %in% genes.intersect)]
# add the disease condition as a new column 
intersect.rfs.final.df.dc <- cbind(intersect.rfs.final.df, condition.disease)

# COMPLETE DATSET HEATMAP #

arr.intersect.rfs.final.df.dc <- intersect.rfs.final.df.dc %>% 
  arrange(factor(rownames(intersect.rfs.final.df.dc), 
              levels = samples.arrange))

# ## melt the dataframe
arr.intersect.rfs.final.df.dc  <- arr.intersect.rfs.final.df.dc %>% select(-condition)
melt.arr.intersect.rfs.final.df.dc <- melt(arr.intersect.rfs.final.df.dc)
melt.arr.intersect.rfs.final.df.dc$samples <- rep(row.names(arr.intersect.rfs.final.df.dc),
                                                  length(arr.intersect.rfs.final.df.dc))



myfilter.intersect <- ifelse(melt.arr.intersect.rfs.final.df.dc$samples %in% nft, "blue", "red")
# blue - normal fallopian tube
# red - omental tumor

pdf("heatmap_intersect.pdf", width = 12, height = 35)
ggplot(melt.arr.intersect.rfs.final.df.dc, aes(samples, variable)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low="white", high="black")  +
  theme(axis.text.x = element_markdown(angle = 90, size = 8, color = myfilter.intersect), 
        axis.text.y = element_text(size=9)) +
  labs(y= "Genes", x="Samples") + theme(legend.position = "none") 
dev.off()


# PLOT SEPARATE HEATMAPS BASED ON DISEASE CONDITION 

# divide the dataframe based on disease condition
rf2.final.nor.c <- rf2.final.df.dc[rf2.final.df.dc$condition == "normal_Fallopian_tube", ]
rf2.final.dis.c <- rf2.final.df.dc[!rf2.final.df.dc$condition== "normal_Fallopian_tube", ]


# arrange the samples 
rf2.final.nor.c <- rf2.final.nor.c %>% arrange(factor(rownames(rf2.final.nor.c), 
                  levels =c('OvCa1850.input', 'OvCa1850.IP', 'OvCa1917.input',
                 'OvCa1917.IP', 'OvCa2005.input', 'OvCa2005.IP',
                 'OvCa2013.input', 'OvCa2013.IP', 
                 'OvCa2053.input', 'OvCa2053.IP', 'OvCa2064.input', 
                 'OvCa2064.IP', 'OvCa2072.input', 'OvCa2072.IP')))

rf2.final.dis.c <- rf2.final.dis.c %>% arrange(factor(rownames(rf2.final.dis.c), 
                  levels =c('OvCa2380.IP','OvCa2380.input', 'OvCa2343.IP',
                          'OvCa2343.input', 'OvCa2270.IP', 'OvCa2270.input',
                          'OvCa2261.IP', 'OvCa2261.input', 'OvCa2221.IP', 
                          'OvCa2221.input','OvCa2186.IP', 'OvCa2186.input')))

# melt the dataframes 
rf2.final.nor.c <- rf2.final.nor.c[,-285]
melt.nor.c <- melt(rf2.final.nor.c)
melt.nor.c$samples <- rep(row.names(rf2.final.nor.c),length(rf2.final.nor.c))


rf2.final.dis.c <- rf2.final.dis.c[,-285]
melt.dis.c <- melt(rf2.final.dis.c)
melt.dis.c$samples <- rep(row.names(rf2.final.dis.c), length(rf2.final.dis.c))


pdf("heatmap_nft.pdf", width = 8, height = 35)
ggplot(melt.nor.c, aes(samples, variable)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 8), axis.text.y = element_text(size=9))+
  labs(y= "Genes", x="Samples") + theme(legend.position = "none") 
dev.off()


pdf("heatmap_ot.pdf", width = 8, height = 35)
ggplot(melt.dis.c, aes(samples, variable)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") + 
  theme(axis.text.x = element_text(angle = 90, size = 8), axis.text.y = element_text(size=9))+
  labs(y= "Genes", x="Samples") + theme(legend.position = "none") 
dev.off()
