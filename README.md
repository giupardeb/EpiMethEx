### HOW-TO
1. scaricare i tre dataset dai seguenti link:

 * SKCM gene expression (pancan normalized) https://genome-cancer.ucsc.edu/download/public/TCGA_SKCM_exp_HiSeqV2_PANCAN-2015-02-24.tgz

 * SKCM DNA methylation Methylation450k https://genome-cancer.ucsc.edu/download/public/TCGA_SKCM_hMethyl450-2015-02-24.tgz 

 * Probset https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL13534&id=11288&db=GeoDb_blob92

2. All'interno del file Script.R settare le seguenti variabili: 

  * `nameFolderDest` sarà il path assoluto dove verranno inseriti i risultati di ogni esecuzione
  
  * `nameFD` cartella che conterrà il range dei geni considerati per ogni esecuzione
  
  * `PathDatasetPancan` path assoluto del dataset **genomicMatrix** all'interno della cartella **pancan_normalized**
  
  * `PathDatasetProbe` path assoluto del dataset **GPL13534-11288.txt** all'interno della cartella **probe**
  
  * `PathDatasetMethylation`path assoluto del dataset **genomicMatrix** all'interno della cartella **Methylation450k**
  
3. Nella riga 30, cambiare il range, all'interno delle parentesi quadre, con quello desiderato. 

**_MAX range di 1000 geni per una congifurazione di PC con 32 GB di RAM e 8 cores_**
