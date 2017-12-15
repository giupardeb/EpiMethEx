## EpiMethEx (Epigenetic Methylation and Expression)

### HOW-TO execute the script

1. Download the following dataset:

 * SKCM gene expression (pancan normalized) https://genome-cancer.ucsc.edu/download/public/TCGA_SKCM_exp_HiSeqV2_PANCAN-2015-02-24.tgz

 * SKCM DNA methylation Methylation450k https://genome-cancer.ucsc.edu/download/public/TCGA_SKCM_hMethyl450-2015-02-24.tgz 

 * Probset https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL13534&id=11288&db=GeoDb_blob92

2. Within the file **Script.R** you must set the following variables: 

  * `nameFolderDest` it's the absolute directory path where will be created within the `nameFD` of every execution.
  
  * `nameFD` it's the directory name which will contain the excel file of every execution.
  
  * `PathDatasetPancan` absolute path of dataset **genomicMatrix** within the directroy **TCGA_SKCM_exp_HiSeqV2_PANCAN-2015-02-24**
  
  * `PathDatasetProbe` absolute path of dataset **GPL13534-11288.txt** within the direcotry **probe**
  
  * `PathDatasetMethylation`path assoluto del dataset **genomicMatrix** all'interno della cartella **TCGA_SKCM_hMethyl450-2015-02-24**
  
3.  Inside the instruction `dfPancan <-fread(PathDatasetPancan,header = T,sep = "\t")[1:10, ]` within square brackets change the range genes

**_MAX range of genes: 1000 with PC configuration: RAM: 32GB CORES: 8_**

### HOW-TO FILTER DATA
