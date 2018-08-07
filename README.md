## EpiMethEx (Epigenetic Methylation and Expression)
The EpiMethEx package is under submission on bioconductor (https://github.com/Bioconductor/Contributions/issues), to use it you can follow the instructions below.
### HOW-TO execute the script
1. clone the repository
2. execute `R CMD build EpiMethEx`
3. install package into R: `install.packages(path_to_file, repos = NULL, type="source")`
4. load library into R: `library(EpiMethEx)`
5. create three datasets as in the example:

```R
 Annotations <- data.frame(
     ID = c("cg1","cg2", "cg3"),
     Relation_to_UCSC_CpG_Island = c("N_Shore","S_Shelf","N_Shore"),
     UCSC_CpG_Islands_Name = c("chrY:11-22","chrY:22-33","chrY:55-66"),
     UCSC_RefGene_Accession = c("NM_00116441;NR_00155","NM_00468","NR_00292"),
     Chromosome_36 = c("Y","Y","Y"),
     Coordinate_36 = c("9973356","21151183","20123684"),
     UCSC_RefGene_Name = c("ARHGEF10L;HIF3A","HIF3A","RNF17"),
     UCSC_RefGene_Group = c("Body;TSS1500","Body","Body"),
     stringsAsFactors=FALSE)
```
```R
 Expressions <- data.frame(
     sample = c("ARHGEF10L", "HIF3A", "RNF17"),
     TCGA1 = c(-0.746592469762, -0.753826336325, -0.531035005853),
     TCGA2 = c(-0.363492469762, -2.67922633632, -0.531035005853),
     TCGA3 = c(-2.97279246976, 0.679073663675, -4.06262633632),
     TCGA4 = c(-0.613189285993,-0.545547600468, 0.59499523215),
     TCGA5 = c(-0.602787155677, -1.56558715568, 3.11398609442),
     TCGA6 = c(0.323339547405, 0.268878821693, -1.62560252829),
     stringsAsFactors=FALSE)
```
```R
 Methylation <- data.frame(
     sample = c("cg1", "cg2", "cg3"),
     TCGA1 = c(0.1290, 0.2394, 0.0913),
     TCGA2 = c(0.9214, 0.0193, 0.0281),
     TCGA3 = c(0.5925, 0.8948, 0.0720),
     TCGA4 = c(0.9643, 0.7934, 0.0354),
     TCGA5 = c(0.1298, 0.0243, 0.8296),
     TCGA6 = c(0.8508, 0.8952, 0.9893),stringsAsFactors=FALSE)
```
5.1 or use the "curatedTCGAData" package:
```R
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite("curatedTCGAData")
library(curatedTCGAData)
library(MultiAssayExperiment)

Methylation <- curatedTCGAData(diseaseCode = "SKCM", assays = "Methylation", dry.run = F)
Expressions <- curatedTCGAData(diseaseCode = "SKCM", assays = "RNASeq2GeneNorm", dry.run = F)
```
it's most important to remember that curatedTCGAData doesn't allow to download dataset of Annotations,therefore it must be loaded manually through csv file or created ad hoc

5.2 or use the csv file:

```R
Expression <- read.csv2("Expressions.csv", header = T,sep = ";",stringsAsFactors=FALSE)
Annotations <- read.csv2("Annotations.csv",header = T,sep = ";",stringsAsFactors=FALSE)
Methylation <- read.csv2("Methylation.csv",header = T,sep = ";",stringsAsFactors=FALSE)
```
6. Execute the following istruction:
 ```R
 epimethex.analysis(Expressions, Annotations, Methylation, 1, 5, 2,TRUE, TRUE, FALSE)
 ```

# Parameters
1. first parameter is a genes expression of data
2. Second parameter is Annotations of data
3. Third parameter is Methylation of data
4. Fourth parameter is a lowerbound of genes
5. fifth parameter is a upperbound of genes
6. sixth parameter, is the number of cores that you can use to analysis
7. seventh parameter determines if genes expression of data are linear
8. eighth parameter determines the test to apply on expression dataset. If TRUE will apply t-student test, otherwise will apply Kolmogorov-Smirnov test
9. ninth parameter, determines the test to apply on methylation dataset. If TRUE will apply t-student test, otherwise will apply Kolmogorov-Smirnov tes

### TESTING (System Configuration)
* OS: Ubuntu 18.04 - 64 bit
* CPU: 8 Cores
* RAM: 32GB
* R: v3.4.4
* MAX range of genes: 1000
