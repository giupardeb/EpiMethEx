#Download Files Tumor

datasetDownload <- function(idMethylation, idPancan){
  
  #create a output folder to files excel
  dir.create(file.path(getwd(), "outputs"), showWarnings = FALSE)
  
  #create folder of datasets
  dir.create(file.path(getwd(), "datasets"), showWarnings = FALSE)
  
  #access into datasets folder
  setwd(file.path(getwd(), "datasets"))
  
  #get name tumor from idMethylation parametr
  tumorName <- unlist(strsplit(idMethylation, ".", fixed = T))[2]
  
  dir.create(file.path(getwd(), tumorName), showWarnings = FALSE)
  
  setwd(file.path(getwd(), tumorName))
  
  nameFileMethylation <- unlist(strsplit(idMethylation, '/', fixed = T))[2]
  nameFilePancan <- unlist(strsplit(idPancan, '/', fixed = T))[2]
  
  #check if not exists gz files mehtylation and not exist methylation file

  checkExistFile(nameFileMethylation, idMethylation)
  
  #check if not exists gz files mehtylation and not exist methylation file
  
  checkExistFile(nameFilePancan, idPancan)
  
  a<-list.files(getwd())
  a[4] <- tumorName
  setwd('..')
  setwd('..')
  return (a)
}

checkExistFile <- function(nameFile, id){
  if(!file_test('-f',paste(nameFile,'.gz',sep = ""))  && !file_test('-f',nameFile)) {
    
    #download file
    baseUrl <- "https://tcga.xenahubs.net/download/"
    
    url <- paste(baseUrl,id,'.gz',sep = "")
    nameFile <- paste(nameFile,'.gz',sep = "")
    download.file(url = url, destfile = nameFile, method = 'wget')
  }
}
