start.time <- Sys.time()

load("/home/giuseppe/github/AnalisiMultigenica/AnalisiMultigenica1giugno.Rdata")
library(jsonlite)
library(mailR)

args <- commandArgs(TRUE)

FCuser <- as.double(args[1])
p_valueUser <- as.double(args[2])

source("functions.R")
#il which ritorna 1,2,3 ciò su dfpancan2 si traducono come: 1=480, 2=481, 3=483

DFfilterGene <- data.frame(dfPancan2[480:482,] >= FCuser)

if (p_valueUser != 0) {
  DFfilterGene1 <-
    data.frame(dfPancan2[c(484, 486, 488),] <= p_valueUser)
  
  #SI POTREBBE SOSTITUIRE ALL'INTERNO DI ANALISI TUTTI GLI NA IN FALSE
  DFfilterGene1[is.na(DFfilterGene1)] <- FALSE
  
  DFfilterGene <- as.data.frame(DFfilterGene1 & DFfilterGene)
  
  DFfilterGene <-
    DFfilterGene[, which(!apply(DFfilterGene == F, 2, all))]
  
} else{
  DFfilterGene <-
    DFfilterGene[, which(!apply(DFfilterGene == F, 2, all))]
}

if (file.exists("geniSelezionti.txt")) {
  GeniSelezionati <- readLines("geniSelezionti.txt")
}


if (exists("GeniSelezionati")) {
  GeniNonSelezionati <-
    GeniSelezionati[which(!GeniSelezionati %in% names(DFfilterGene))]
  
  Geni <-
    GeniSelezionati[which(GeniSelezionati %in% names(DFfilterGene))]
  
  #filtro il dataframe in base ai geni selezionati dall'utente
  DFfilterGene <-
    as.data.frame(DFfilterGene[, Geni])
  colnames(DFfilterGene) <- Geni
  
}


dimDF <- ncol(DFfilterGene)
outputJson <- list()

if (dimDF != 0) {
  for (i in 1:dimDF) {
    #aggiungo un gene nel JSON
    numGene <<- paste("gene", i, sep = "")
    outputJson[[numGene]][["CG"]] <- list()
    outputJson[[numGene]] <- list()
    
    
    namesGene <- colnames(DFfilterGene[, i, drop = F])
    
    if (DFfilterGene[1, i] == "TRUE" &&
        DFfilterGene[2, i] == "TRUE" &&
        DFfilterGene[3, i] == "TRUE") {
      #ELABORO TUTTO
      
      getValueDFPancan2(namesGene, c(477:488))
      getValueDFFinale(namesGene, c(1:6))
      
    }
    else if (DFfilterGene[1, i] == "TRUE" &&
             DFfilterGene[2, i] == "TRUE") {
      #ELABORO UPvsMID E UPvsDOWN
      
      getValueDFPancan2(namesGene, c(477:481, 483:486))
      getValueDFFinale(namesGene, c(1:4))
      
    }
    else if (DFfilterGene[1, i] == "TRUE" &&
             DFfilterGene[3, i] == "TRUE") {
      #ELABORO UPvsMID E MIDvsDOWN
      
      getValueDFPancan2(namesGene, c(477:480, 482:484, 487, 488))
      getValueDFFinale(namesGene, c(1:2, 5:6))
    }
    else if (DFfilterGene[2, i] == "TRUE" &&
             DFfilterGene[3, i] == "TRUE") {
      #ELABORO UPvsDOWN E MIDvsDOWN
      
      getValueDFPancan2(namesGene, c(477:479, 481, 482, 485:488))
      getValueDFFinale(namesGene, c(3:6))
    }
    else if (DFfilterGene[1, i] == "TRUE") {
      #ELABORO UPvsMID
      
      getValueDFPancan2(namesGene, c(477:480, 483, 484))
      getValueDFFinale(namesGene, c(1:2))
    }
    else if (DFfilterGene[2, i] == "TRUE") {
      #ELABORO UPvsDOWN
      
      getValueDFPancan2(namesGene, c(477:479, 481, 485, 486))
      getValueDFFinale(namesGene, c(3:4))
    }
    else if (DFfilterGene[3, i] == "TRUE") {
      #ELABORO MIDvsDOWN
      
      getValueDFPancan2(namesGene, c(477:479, 482, 487, 488))
      getValueDFFinale(namesGene, c(5:6))
    }
  }
}else{
  outputJson['errore'] <- "Il Gene che hai selezionato non soddisfa i parametri da te indicati"
}

#print(toJSON(outputJson, pretty = TRUE, auto_unbox = TRUE))
write(toJSON(outputJson, pretty = TRUE, auto_unbox = TRUE), file = "test.json")

send.mail(
  from = "relief.portal@gmail.com",
  to = "giupardeb@yahoo.it",
  subject = "il tuo Json è pronto!",
  body = "ricevi questa mail poiché il tuo Json è pronto, visualizzalo qui: http://localhost/result.php",
  smtp = list(
    host.name = "smtp.gmail.com",
    port = 465,
    user.name = "relief.portal",
    passwd = "",
    ssl = TRUE
  ),
  authenticate = TRUE,
  send = TRUE
)
end.time <- Sys.time()
print(end.time - start.time)
