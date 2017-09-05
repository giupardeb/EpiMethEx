library(readxl)
library(data.table)
library(ggplot2)
require(reshape2)
require(ggplot2)

A <-
  read_excel("/home/giuseppe/github/AnalisiMultigenica/dataset/analisiMMP9.xlsx")

#boxplot cg global
valori <- A[14:486, 4:23]
valori <- apply(valori, 2, as.double)
valori.m <- melt(valori)
valori.m$Var2 <- factor(valori.m$Var2, unique(as.character(valori.m$Var2)))
ggplot(data=valori.m,aes(x=Var2,y=value))+ geom_boxplot(outlier.shape = NA,fill='#619cff', color="black")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot MMP9
valori <- A[14:486,2:3]
valori.m <- melt(valori, id.vars = "X__2")
valori.m$value<-as.numeric(as.character(valori.m$value))
valori.m$X__2 <- factor(valori.m$X__2, unique(as.character(valori.m$X__2)))

ggplot(data = valori.m, aes(x = X__2, y = value)) + geom_boxplot(aes(fill = X__2),outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_manual(name = "MMP9", values = c("#619cff", "#00ba38","#f8766d"))

#boxplot cg con up mid down uniti
valori <- A[14:486, 2:23]
valori <- subset(valori, select = -MMP9)
valori.m <- melt(valori, id.vars = "X__2")
valori.m$value<--as.numeric(as.character(valori.m$value))
valori.m$X__2 <- factor(valori.m$X__2, unique(as.character(valori.m$X__2)))
ggplot(data = valori.m, aes(x = variable, y = value)) + geom_boxplot(aes(fill = X__2),outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_manual(name = "MMP9", values = c("#619cff", "#00ba38","#f8766d"))

#boxplot CG per posizione

valori <- A[11:486,4:23]
valori<-valori[-2:-3,]
df<-matrix(NA,nrow = 9460,ncol = 3)
colnames(df)<-c('cg','posizione','value')
posizioni<-as.character(valori[1,])
cg<-colnames(valori)
valori<-t(valori)
valori <- valori[,-1]

i<-1
  for(j in 1:20){#identifica cg e posizioni
    for(k in 1:473){
      df[i,1]<-cg[j]
      df[i,2]<-posizioni[j]
      df[i,3]<-valori[j,k]
      i<-i+1
    }
  }
df<-as.data.frame(df)
df$value<-as.numeric(as.character(df$value))
df<-cbind(df,A[14:486,2])
df$f1f2 <- interaction(df$cg, df$X__2)

#df$f1f2 <- factor(df$f1f2, unique(as.character(df$f1f2)))
df$X__2 <- factor(df$X__2, unique(as.character(df$X__2))) #serve per avere l'ordine UP MID DOWN
df$posizione <- factor(df$posizione, unique(as.character(df$posizione)))
ggplot(data = df, aes(x =posizione , y = value)) + geom_boxplot(aes(fill = factor(X__2)),outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_manual(name = "MMP9", values = c("#619cff", "#00ba38","#f8766d",'#b978a3'))

#boxplot CpG_Island_pos
valori <- A[12:486,4:23]
valori<-valori[-2,]
AA<-t(valori)
AA<-as.data.frame(AA)

s <- strsplit(as.character(AA$V1), split = ":") #da inserire
GPL2 <-
  data.frame(V1 = rep(colnames(valori), sapply(s, length)), V22 = unlist(s))#da inserire

nrows <- 473 * 40 #473=valori cg, 40 numero di righe di GPL2

df<-matrix(NA,nrow = nrows,ncol = 3)#da inserire
colnames(df)<-c('cg','identificativo','value')#da inserire
posizioni<-as.character(GPL2$V22)#da inserire
cg<-colnames(valori)#da inserire
valori<-t(valori)
valori <- valori[,-1]

z<-1
i<-1
for(j in rep(1:20, each = 2)){
  for(k in 1:473){
    df[i,1]<-cg[j]
    df[i,2]<-posizioni[z]
    df[i,3]<-as.numeric(valori[j,k])
    i<-i+1
  }
  z<-z+1
}
df<-as.data.frame(df)
#df$value<-as.numeric(df$value)
df<-cbind(df,A[14:486,2])
df$f1f2 <- interaction(df$cg, df$X__2)

df$value<-as.numeric(as.character(df$value))
df$identificativo <- factor(df$identificativo, unique(as.character(df$identificativo)))
df$X__2 <- factor(df$X__2, unique(as.character(df$X__2)))
ggplot(data = df, aes(x = identificativo, y = value)) + geom_boxplot(aes(fill = factor(X__2)),outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_fill_manual(name = "MMP9", values = c("#619cff", "#00ba38","#f8766d",'#b978a3'))
