# Crosstalk between miRNA expression and DNA methylation drives the hormone-dependent phenotype of breast cancer

*Aure MR, Fleischer T, Bjørklund S, Ankill J, Castro-Mondragon JA; OSBREAC, Børresen-Dale AL, Tost J, Sahlberg KK, Mathelier A, Tekpli X, Kristensen VN.*
Genome Med. 2021 Apr 29;13(1):72. doi: 10.1186/s13073-021-00880-4.
PMID: 33926515 

Contact: Miriam Ragle Aure,  m.r.aure@medisin.uio.no

The code and associated data is presented as a minimal working example for miRNA-methylation Quantitative Trait Loci (mimQTL) analysis presented in the manuscript referenced above.

In brief, we take two matrices with data from the same samples, mirMat holding normalized miRNA expression data, and metMat holding normalized CpG methylation data, and correlate all rows in the matrices (e.g. all miRNAs to all CpGs). The significant miRNA-CpG associations, i.e. mimQTLs, are kept and clustered. The final number of miRNA and CpG clusters needs to be decided by manual inspection. Here, as an example, two miRNA clusters and two CpG clusters are chosen. These and default parameters such as cut-off points for filtering CpGs with low variation across samples or miRNAs expressed in few samples should be adapted to the data.

The present example consists of 100 random samples with 1000 random CpGs and 100 random miRNAs taken from The Cancer Genome Atlas breast cancer (TCGA-BRCA) data (The Cancer Genome Atlas Network: Comprehensive molecular portraits of human breast tumours. Nature 2012, 490:61-70) downloaded from the Xena browser on April 6th 2020 (Goldman M et al.:The UCSC Xena platform for public and private cancer genomics data visualization and interpretation. bioRxiv 2019:326470), https://xenabrowser.net/datapages/?cohort=TCGA%20Breast%20Cancer%20(BRCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
