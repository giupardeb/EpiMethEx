#Download tumor files

datasetDownload <- function(idMethylation, idPancan){
  
  #create output folder for excel files
  dir.create(file.path(getwd(), "outputs"), showWarnings = FALSE)
  
  #create folder for datasets
  dir.create(file.path(getwd(), "datasets"), showWarnings = FALSE)
  
  #access  datasets folder
  setwd(file.path(getwd(), "datasets"))
  
  #get tumor name from idMethylation parameter
  tumorName <- unlist(strsplit(idMethylation, ".", fixed = T))[2]
  
  #create a folder with the name of the tumor
  dir.create(file.path(getwd(), tumorName), showWarnings = FALSE)
  
  #access "tumorName" folder
  setwd(file.path(getwd(), tumorName))
  
  nameFileMethylation <- unlist(strsplit(idMethylation, '/', fixed = T))[2]
  nameFilePancan <- unlist(strsplit(idPancan, '/', fixed = T))[2]
  
  #check the non existence of gz files mehtylation and  methylation file
  checkExistFile(nameFileMethylation, idMethylation)
  
  #check the non existence of gz files mehtylation and methylation file
  checkExistFile(nameFilePancan, idPancan)
  
  #the variable "params" will contain the following files:
  # params[1] = "methylation.gz"
  # params[2] = "pancan.gz"
  # params[3] = "probeset.txt"
  # params[4] = "tumorName"
  params <- list.files(getwd())
  
  #add "tumorName" into params
  params[4] <- tumorName
  
  setwd('..')
  setwd('..')
  
  return (params)
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
