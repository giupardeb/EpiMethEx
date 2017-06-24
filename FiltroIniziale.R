FCuser <- 2
p_valueUser <- 0.05
limit <- 100
source("functions.R")
library(jsonlite)
library(dplyr)
load("/home/giuseppe/github/AnalisiMultigenica/AnalisiMultigenica1giugno.Rdata")

UPvsMIDordinato <- dfPancan2[c("fc_UPvsMID", "pvalue_UPvsMID"), ]
UPvsDOWNordinato <- dfPancan2[c("fc_UPvsDOWN", "pvalue_UPvsDOWN"), ]
MIDvsDOWNordinato <-
  dfPancan2[c("fc_MIDvsDOWN", "pvalue_MIDvsDOWN"), ]

UPvsMIDordinato <-
  UPvsMIDordinato[c(as.numeric(which(UPvsMIDordinato[1,] >= FCuser)))]
UPvsDOWNordinato <-
  UPvsDOWNordinato[c(as.numeric(which(UPvsDOWNordinato[1,] >= FCuser)))]
MIDvsDOWNordinato <-
  MIDvsDOWNordinato[c(as.numeric(which(MIDvsDOWNordinato[1,] >= FCuser)))]


UPvsMIDordinato <-
  UPvsMIDordinato[c(as.numeric(which(UPvsMIDordinato[2,] <= p_valueUser)))]
UPvsDOWNordinato <-
  UPvsDOWNordinato[c(as.numeric(which(UPvsDOWNordinato[2,] <= p_valueUser)))]
MIDvsDOWNordinato <-
  MIDvsDOWNordinato[c(as.numeric(which(MIDvsDOWNordinato[2,] <= p_valueUser)))]


# #raggruppa i cg ma non raggruppa FC GENE PVALUE
# Z.Z<-C %>%
#   group_by(Gene=C$Gene,FC=C$FC,pValue=C$pValue,CG=C$CG) %>%
#   summarise(Posizione=paste(Identificativo,Posizione,collapse=";;"))

# ordinamento(UPvsMIDordinato, "UPvsMIDordinato")
# ordinamento(UPvsDOWNordinato, "UPvsDOWNordinato")
# ordinamento(MIDvsDOWNordinato, "MIDvsDOWNordinato")



A <- as.data.frame(t(UPvsDOWNordinato))
colnames(A) <- c("FC", "pValue")
A <- cbind(A, V22 = rownames(A))

if (!(nrow(A) < limit)) {
  
  A <- A[order(A$FC, decreasing = T), ]
  A <- rbind(head(A, limit / 2), tail(A, limit / 2))
  
}

C <- merge(A, righeCheTiServono1, sort = F)

colnames(C) <-
  c("Gene", "FC", "pValue", "CG", "Identificativo", "Posizione")

#raggruppa FC GENE PVALUE, i cg sono all'interno di posizione
Z<-C %>%
  group_by(Gene=C$Gene,FC=C$FC,pValue=C$pValue) %>%
  summarise(Posizione=paste(CG,Identificativo,Posizione,sep=";;",collapse=";;"))

UPvsDOWNordinato <- Z[order(Z$FC, decreasing = T), ]


A <- as.data.frame(t(UPvsMIDordinato))
colnames(A) <- c("FC", "pValue")
A <- cbind(A, V22 = rownames(A))

if (!(nrow(A) < limit)) {
  
  A <- A[order(A$FC, decreasing = T), ]
  A <- rbind(head(A, limit / 2), tail(A, limit / 2))
  
}

C <- merge(A, righeCheTiServono1, sort = F)

colnames(C) <-
  c("Gene", "FC", "pValue", "CG", "Identificativo", "Posizione")

#raggruppa FC GENE PVALUE, i cg sono all'interno di posizione
Z<-C %>%
  group_by(Gene=C$Gene,FC=C$FC,pValue=C$pValue) %>%
  summarise(Posizione=paste(CG,Identificativo,Posizione,sep=";;",collapse=";;"))

UPvsMIDordinato <- Z[order(Z$FC, decreasing = T), ]



A <- as.data.frame(t(MIDvsDOWNordinato))
colnames(A) <- c("FC", "pValue")
A <- cbind(A, V22 = rownames(A))

if (!(nrow(A) < limit)) {
  
  A <- A[order(A$FC, decreasing = T), ]
  A <- rbind(head(A, limit / 2), tail(A, limit / 2))
  
}

C <- merge(A, righeCheTiServono1, sort = F)

colnames(C) <-
  c("Gene", "FC", "pValue", "CG", "Identificativo", "Posizione")

#raggruppa FC GENE PVALUE, i cg sono all'interno di posizione
Z<-C %>%
  group_by(Gene=C$Gene,FC=C$FC,pValue=C$pValue) %>%
  summarise(Posizione=paste(CG,Identificativo,Posizione,sep=";;",collapse=";;"))

MIDvsDOWNordinato <- Z[order(Z$FC, decreasing = T), ]


write(toJSON(unname(split(
  UPvsMIDordinato, 1:nrow(UPvsMIDordinato)
))), file = "UPvsMIDordinato.json")

write(toJSON(unname(split(
  UPvsDOWNordinato, 1:nrow(UPvsDOWNordinato)
))), file = "UPvsDOWNordinato.json")

write(toJSON(unname(split(
  MIDvsDOWNordinato, 1:nrow(MIDvsDOWNordinato)
))), file = "MIDvsDOWNordinato.json")
