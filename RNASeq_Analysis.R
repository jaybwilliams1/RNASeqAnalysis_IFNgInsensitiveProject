library(limma)
library(edgeR)
library(stringr)

setwd("~/Documents/alignments")

#Import files
i1 <- read.table("IFNgR2KO_Mu1_abundance.tsv")
i2 <- read.table("IFNgR2KO_Mu2_abundance.tsv")
i3 <- read.table("IFNgR2KO_Mu3_abundance.tsv")

j1 <- read.table("Jak1KO_Mu1_abundance.tsv")
j2 <- read.table("Jak1KO_Mu2_abundance.tsv")
j3 <- read.table("Jak1KO_Mu3_abundance.tsv")

w1 <- read.table("WT_Mu1_abundance.tsv")
w2 <- read.table("WT_Mu2_abundance.tsv")
w3 <- read.table("WT_Mu3_abundance.tsv")

#Select columns with target_id and est_counts
i1 <- i1[,c("V1","V4")]
i2 <- i2[,c("V1","V4")]
i3 <- i3[,c("V1","V4")]

j1 <- j1[,c("V1","V4")]
j2 <- j2[,c("V1","V4")]
j3 <- j3[,c("V1","V4")]

w1 <- w1[,c("V1","V4")]
w2 <- w2[,c("V1","V4")]
w3 <- w3[,c("V1","V4")]

#Rename counts column to match sample
colnames(i1) <- c("gene_id", "i1")
colnames(i2) <- c("gene_id", "i2")
colnames(i3) <- c("gene_id", "i3")

colnames(j1) <- c("gene_id", "j1")
colnames(j2) <- c("gene_id", "j2")
colnames(j3) <- c("gene_id", "j3")

colnames(w1) <- c("gene_id", "w1")
colnames(w2) <- c("gene_id", "w2")
colnames(w3) <- c("gene_id", "w3")

#Merge samples into one data frame
m1 <- merge(i1, i2, by.x = 1, by.y = 1)
m2 <- merge(m1, i3, by.x = 1, by.y = 1)
m3 <- merge(m2, j1, by.x = 1, by.y = 1)
m4 <- merge(m3, j2, by.x = 1, by.y = 1)
m5 <- merge(m4, j3, by.x = 1, by.y = 1)
m6 <- merge(m5, w1, by.x = 1, by.y = 1)
m7 <- merge(m6, w2, by.x = 1, by.y = 1)
m8 <- merge(m7, w3, by.x = 1, by.y = 1)

#Divide target_id column into parts
FOO <- data.frame(do.call('rbind', strsplit(as.character(m8$gene_id), '|', fixed=TRUE)))

#Select gene_id, gene_name, and gene_type
sub <- FOO[,c("X2", "X6", "X8")]

#Merge selected identifier columns with count data
new <- cbind(sub, m8)

#Subset protein coding genes
protein_coding <- new[new$X8 == "protein_coding",]
row.names(protein_coding) <- protein_coding$gene_id

#Subset gene_id and gene_name
gene_names <- protein_coding[, c("X2", "X6")]

#Select just counts for protein coding genes and save matrix
merged2 <- protein_coding[, c("i1", "i2", "i3", "j1", "j2", "j3", "w1", "w2", "w3")]
write.csv(merged2, "jaycounts.csv")

exp <- read.csv("jaycounts.csv", row.names = 1, fill = T)

#Make DGElist
lcpm <- cpm(exp, log = T)

x <- DGEList(counts=exp)
x$genes <- gene_names

sub <- c("i1", "i2", "i3", "j1", "j2", "j3", "w1", "w2", "w3")
set <- c("1", "1", "1", "2", "2", "2", "3", "3", "3")
subset <- data.frame(sub, set)
x$samples <- merge(x$samples, subset, by.x = 0, by.y = 1, all.x = TRUE)
row.names(x$samples) <- x$samples[,1]
x$samples <- x$samples[,-1]
x$samples <- x$samples[,-1]
names <- c("lib.size", "norm.factors", "group")
colnames(x$samples) <- names
x$samples

#Filter lowly expressed genes
keep.exprs <- rowSums(lcpm>=1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#Normalize
x <- calcNormFactors(x, method = "TMM")

#Create design matrix
design <- model.matrix(~0+group, data=x$samples) 
design
 
#Voom
v=voom(x,design,plot=T)

#Create contrast matrix
contr.matrix <- makeContrasts(
  ivsw = group1 - group3,
  jvsw = group2 - group3,
  levels = colnames(design))
contr.matrix

#Fit model
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
dim(efit)
plotSA(efit, main="Final model: Mean−variance trend")

summary(decideTests(efit))

#Restrict model based on log fold change
tfit  <- treat(vfit, lfc = 1)
dt <- decideTests(tfit)
summary(dt)

#Examine p-value tables
table <- topTable(tfit,coef="ivsw",sort.by="p", n=Inf)
table <- table[,-1]
tablej <- topTable(tfit,coef="jvsw",sort.by="p", n=Inf)
tablej <- tablej[,-1]

write.csv(table, "ivswdeg.csv")
write.csv(tablej, "jvswdeg.csv")

#Make venn diagram of overlapping genes
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n=20)
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))