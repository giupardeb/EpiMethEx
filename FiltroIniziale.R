source("functionsFiltriIniziali.R")
library(jsonlite)
library(dplyr)
load("/home/giuseppe/github/AnalisiMultigenica/AnalisiMultigenica1giugno.Rdata")

FCuser <- 2
p_valueUser <- 0.05
n <- 50

UPvsMIDordinato <- dfPancan2[c("fc_UPvsMID", "pvalue_UPvsMID"),]
UPvsDOWNordinato <- dfPancan2[c("fc_UPvsDOWN", "pvalue_UPvsDOWN"),]
MIDvsDOWNordinato <-
  dfPancan2[c("fc_MIDvsDOWN", "pvalue_MIDvsDOWN"),]

UPvsMIDordinato <-
  UPvsMIDordinato[c(as.numeric(which(
    UPvsMIDordinato[1, ] >= FCuser |
      UPvsMIDordinato[1, ] <= (-FCuser)
  )))]
UPvsDOWNordinato <-
  UPvsDOWNordinato[c(as.numeric(which(
    UPvsDOWNordinato[1, ] >= FCuser |
      UPvsDOWNordinato[1, ] <= (-FCuser)
  )))]
MIDvsDOWNordinato <-
  MIDvsDOWNordinato[c(as.numeric(which(
    MIDvsDOWNordinato[1, ] >= FCuser |
      MIDvsDOWNordinato[1, ] <= (-FCuser)
  )))]


UPvsMIDordinato <-
  UPvsMIDordinato[c(as.numeric(which(UPvsMIDordinato[2, ] <= p_valueUser)))]
UPvsDOWNordinato <-
  UPvsDOWNordinato[c(as.numeric(which(UPvsDOWNordinato[2, ] <= p_valueUser)))]
MIDvsDOWNordinato <-
  MIDvsDOWNordinato[c(as.numeric(which(MIDvsDOWNordinato[2, ] <= p_valueUser)))]

ordinamento(UPvsMIDordinato, "UPvsMIDordinato")
ordinamento(UPvsDOWNordinato, "UPvsDOWNordinato")
ordinamento(MIDvsDOWNordinato, "MIDvsDOWNordinato")


write(toJSON(unname(split(
  UPvsMIDordinato, 1:nrow(UPvsMIDordinato)
))), file = "UPvsMIDordinato.json")

write(toJSON(unname(split(
  UPvsDOWNordinato, 1:nrow(UPvsDOWNordinato)
))), file = "UPvsDOWNordinato.json")

write(toJSON(unname(split(
  MIDvsDOWNordinato, 1:nrow(MIDvsDOWNordinato)
))), file = "MIDvsDOWNordinato.json")

# #raggruppa i cg ma non raggruppa FC GENE PVALUE
# Z.Z<-C %>%
#   group_by(Gene=C$Gene,FC=C$FC,pValue=C$pValue,CG=C$CG) %>%
#   summarise(Posizione=paste(Identificativo,Posizione,collapse=";;"))

