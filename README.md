## EpiMethEx (Epigenetic Methylation and Expression)
The EpiMethEx package is under submission on bioconductor (https://github.com/Bioconductor/Contributions/issues), to use it you can follow the instructions below.
### HOW-TO execute the script
1. clone the repository
2. execute `R CMD build EpiMethEx`
3. install package into R: `install.packages(path_to_file, repos = NULL, type="source")`
4. create three datasets as in the example. 
The example below is a extract of the datasets TCGA-SKCM 
(https://xenabrowser.net/datapages/?cohort=TCGA%20Melanoma%20(SKCM)):

```R
 Annotations <- data.frame(
 ID = c("cg11663302","cg01552731", "cg09081385"),
 Relation_to_UCSC_CpG_Island = c("Island","N_Shore","N_Shore"),
 UCSC_CpG_Islands_Name = c("chr1:18023481-18023792","chr19:46806998-46807617",
 "chr12:120972167-120972447"),
 UCSC_RefGene_Accession = c("NM_001011722","NM_152794","NM_014868"),
 Chromosome_36 = c("1","19","12"),
 Coordinate_36 = c("17896255","51498747","119456453"),
 UCSC_RefGene_Name = c("ARHGEF10L","HIF3A","RNF10"),
 UCSC_RefGene_Group =c("Body","1stExon","TSS200"),
 stringsAsFactors=FALSE)
```
```R
 Expressions <- data.frame(
 'sample' = c("ARHGEF10L", "HIF3A", "RNF10"),
 'TCGA-YD-A89C-06' = c(-0.746592469762, -0.753826336325, 0.4953280),
 'TCGA-Z2-AA3V-06' = c(0.578807530238, -2.30662633632, 0.1023280),
 'TCGA-EB-A3Y6-01' = c(-0.363492469762, -2.67922633632, -0.6147720),
 'TCGA-EE-A3JA-06' = c(-2.97279246976, -3.61932633632, 0.02932801),
 'TCGA-D9-A4Z2-01' = c(-0.128492469762, 0.679073663675, 0.4017280),
 'TCGA-D3-A51G-06' = c(-0.4299925, -4.0626263, -1.0136720),
 stringsAsFactors=FALSE)
```
```R
 Methylation <- data.frame(
 'sample' = c("cg11663302", "cg01552731", "cg09081385"),
 'TCGA-YD-A89C-06' = c(0.9856, 0.7681, 0.0407),
 'TCGA-Z2-AA3V-06' = c(0.9863, 0.8551, 0.0244),
 'TCGA-EB-A3Y6-01' = c(0.9876, 0.6473, 0.028),
 'TCGA-EE-A3JA-06' = c(0.9826, 0.4587, 0.0343),
 'TCGA-D9-A4Z2-01' = c(0.9881, 0.8509, 0.0215),
 'TCGA-D3-A51G-06' = c(0.9774, 0.813, 0.0332),
 stringsAsFactors=FALSE)
```
4.1 or use the "curatedTCGAData" package:
```R
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("curatedTCGAData")
library(curatedTCGAData)
library(MultiAssayExperiment)

Methylation <- curatedTCGAData(diseaseCode = "SKCM", assays = "Methylation", dry.run = F)
Expressions <- curatedTCGAData(diseaseCode = "SKCM", assays = "RNASeq2GeneNorm", dry.run = F)
```
it's most important to remember that curatedTCGAData doesn't allow to download dataset of Annotations,therefore it must be loaded manually through csv file or created ad hoc

4.2 or use the csv file:

```R
Expression <- read.csv2("Expressions.csv", header = T,sep = ";",stringsAsFactors=FALSE)
Annotations <- read.csv2("Annotations.csv",header = T,sep = ";",stringsAsFactors=FALSE)
Methylation <- read.csv2("Methylation.csv",header = T,sep = ";",stringsAsFactors=FALSE)
```
5. Execute the following istruction:
 ```R
 library(EpiMethEx)
 epimethex.analysis(Expressions, Annotations, Methylation, 1, 3, 2,TRUE, TRUE, FALSE)
 ```

# Parameters
1. first parameter is a genes expression of data
2. Second parameter is Annotations of data
3. Third parameter is Methylation of data
4. Fourth parameter and  fifth parameter are the range of genes to consider
5. sixth parameter, is the number of cores that you can use to analysis
6. seventh parameter determines if genes expression of data are linear
7. eighth parameter determines the test to apply on expression dataset. If TRUE will apply t-student test, otherwise will apply Kolmogorov-Smirnov test
8. ninth parameter, determines the test to apply on methylation dataset. If TRUE will apply t-student test, otherwise will apply Kolmogorov-Smirnov test

### HOW TO FILTER

To further evaluate the biological significance of the methylation hotspots involved in gene regulation mechanisms, is possible to filter EpiMethEx output with an additional R script (https://github.com/giupardeb/EpiMethEx-Filter)

### TESTING (System Configuration)
* OS: Ubuntu 18.04 - 64 bit
* CPU: 8 Cores
* RAM: 32GB
* R: v3.4.4
* MAX range of genes: 1000
