###############################################################################
# The code and associated data is presented as a minimal working example for
# miRNA-methylation Quantitative Trait Loci (mimQTL) analysis presented in:

# Title:
# Crosstalk between miRNA expression and DNA methylation drive the hormone-
# dependent phenotype of breast cancer. bioRxiv 2020.

# Authors:
# Miriam Ragle Aure, Thomas Fleischer, Sunniva Bjørklund, Jørgen Ankill,
# Jaime A Castro-Mondragon, OSBREAC, Anne-Lise Børresen-Dale, Kristine
# K. Sahlberg, Anthony Mathelier, Xavier Tekpli and Vessela N. Kristensen.

# Contact:
# Miriam Ragle Aure (2020-April)
# m.r.aure@medisin.uio.no

# In brief, we take two matrices with data from the same samples, mirMat holding
# normalized miRNA expression data, and metMat holding normalized CpG
# methylation data, and correlate all rows in the matrices (e.g. all miRNAs to
# all CpGs). The significant miRNA-CpG associations, i.e. mimQTLs, are kept and
# clustered. The final number of miRNA and CpG clusters needs to be decided by
# manual inspection. Here, as an example, two miRNA clusters and two CpG
# clusters are chosen. These and default parameters such as cut-off points for
# filtering CpGs with low variation across samples or miRNAs expressed in few
# samples should be adapted to the data.

# The present example consists of 100 random samples and 1000 random CpGs and 100 random miRNAs
# of The Cancer Genome Atlas breast cancer (TCGA-BRCA) data downloaded from the
# Xena browser April 6th 2020, (Goldman M et al.:The UCSC Xena platform for
# public and private cancer genomics data visualization and interpretation.
# bioRxiv 2019:326470), https://xenabrowser.net/datapages/?cohort=TCGA%20Breast%20Cancer%20(BRCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

###############################################################################
## dependencies
###############################################################################

library(reshape2) # CRAN
library(pheatmap) # CRAN

###############################################################################
## load data, match column names and check same order
###############################################################################

# miRNA expression
mirMat = data.matrix(read.csv("miRNA_HiSeq_gene_TCGA_100x100.txt",
                    sep="\t", row.names=1))
head(mirMat[,1:3])

# CpG methylation
metMat = data.matrix(read.csv("HumanMethylation450_TCGA_1000x100.txt",
                              sep="\t", row.names=1))
head(metMat[,1:3])

sampleNames = intersect(colnames(mirMat), colnames(metMat))
mirMat = mirMat[,sampleNames]
metMat = metMat[,sampleNames]
all(colnames(metMat)==colnames(mirMat))

dim(mirMat) # 100 x 100
dim(metMat) # 1000 x 100

###############################################################################
## feature filtering
###############################################################################

# filter out miRNAs that are detected in <=10% of the samples
m = min(mirMat)
det = as.matrix(mirMat[,] > m)
dm = apply(det, 1, mean)
mirMat = mirMat[dm>0.1, ]

# filter out CpGs with low variation
iqr = apply(metMat, 1, IQR)
metMat = metMat[iqr>0.1,]

###############################################################################
## see data distributions
###############################################################################

par(mfrow=c(1,2))
hist(metMat, main="methylation_TCGA_subset")
hist(mirMat, main="mirna_TCGA_subset")

###############################################################################
## correlate miRNA and methylation with Spearman correlation
###############################################################################

tab = matrix(data=NA, nrow=nrow(metMat), ncol = nrow(mirMat))
pvaltab = matrix(data=NA, nrow=nrow(metMat), ncol = nrow(mirMat))

# ~2 minutes run-time with present data set on a 2016 laptop
# expect warnings for exact p-values and ties
for (i in 1:nrow(metMat)){
  for (j in 1:nrow(mirMat)) {
    testres = cor.test(metMat[i, ], mirMat[j, ], method="spearman")
    tab[i,j] = testres$estimate
    pvaltab[i,j] = testres$p.value
  }
}

rownames(tab) = rownames(metMat)
colnames(tab) = rownames(mirMat)
rownames(pvaltab) = rownames(metMat)
colnames(pvaltab) = rownames(mirMat)

###############################################################################
## keep significant correlations after Bonferroni correction
###############################################################################

par(mfrow=c(1,1))

# expect non-uniform distribution
hist(pvaltab, xlab = "correlation p-values")

keep = pvaltab < 0.05/(nrow(pvaltab)*ncol(pvaltab))
table(keep)

sig.tab = c()
for (i in 1:nrow(pvaltab)) {
  for (j in 1:ncol(pvaltab)) {
    if (pvaltab[i,j] < 0.05/(nrow(pvaltab)*ncol(pvaltab))) {
      sig.tab = append(sig.tab,
             list(data.frame(
                CpG=rownames(pvaltab)[i],
                miRNA=colnames(pvaltab)[j],
                Spearman_pval=pvaltab[i,j],
                Spearman_cor=tab[i,j])))
    }
  }
}

sig.tab <- do.call(rbind, sig.tab)
head(sig.tab)
dim(sig.tab)

###############################################################################
## making a matrix of significant p-values with pos or neg correlation indicated
###############################################################################

# make a matrix object out of input data, e.g. all significant CpGs in rows and
# all significant miRNAs in columns
d = dcast(sig.tab, CpG ~ miRNA, value.var = "Spearman_cor")
m = as.matrix(d[,-1])
rownames(m) = d[,1]

# convert to a -1/0/1 matrix
m[m>0] =  1
m[m<0] = -1
m[is.na(m)] = 0

###############################################################################
## Clustering with pheatmap
###############################################################################

mycolors = colorRampPalette(c("dodgerblue3","white", "firebrick1"))(n=299)

# Correlation distance and average linkage
pheatmap(m, clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         cutree_rows=2, cutree_cols=2, # decide number of clusters
         annotation_colors = annotations_colors,
         col = mycolors,
         show_rownames=TRUE,
         show_colnames=TRUE)

###############################################################################
## get and save pheatmap clusters
###############################################################################

out = pheatmap(m,
        clustering_distance_rows = "correlation",
        clustering_distance_cols = "correlation",
        clustering_method = "average")

# CpGs in rows
a = rownames(m[out$tree_row[["order"]],])
b = sort(cutree(out$tree_row, k=2))
write.table(b, "CpG_clusters.txt", sep = "\t")

# miRNAs in columns
c = rownames(m[out$tree_col[["order"]],])
d = sort(cutree(out$tree_col, k=2))
write.table(d, "miRNA_clusters.txt", sep = "\t")

###############################################################################
sessionInfo()
