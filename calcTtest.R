#args <- commandArgs(TRUE) args dovr? contenere il path inviato da php
start.time <- Sys.time()
library(RJSONIO)
library(data.table)

#questi dataframe servono per la seconda parte quando si va alla ricerca delle etichette
df = fread("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/probe/GPL13534-11288.txt", sep="\t",header=F,skip=37)

df1 = fread("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/Methylation450k/genomicMatrix", sep="\t")
#df2 = fread("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/pancan_normalized/genomicMatrix", sep="\t")
conn <- file("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/pancan_normalized/genomicMatrix", open="r") 
line <-readLines(conn)
close(conn)

errori = 0
data = I
json1 <- list()



dataIsLinear = T #da chiedere all'utente
parametro = 0.05 #utilizzato per il foldchange, da chiedere all'utente
etichette = tail(strsplit(line[1],"\t")[[1]], -1)

#vado alla ricerca delle etichette nel dataset GPL e 
#succ. calcolerò i ttest trovati nel terzo dataset
righeCheTiServono<- df[which(df$V22=="MMP9"), ]



dataFrameFinale = NULL
for(i in 1:length(righeCheTiServono$V1)){
  
  dataFrame<- df1[which(df1$sample==righeCheTiServono$V1[i]), ]
  dataFrameFinale = rbind(dataFrameFinale,dataFrame[1, ])
}

dataFrameFinale = na.omit(dataFrameFinale)



#length(line)
for(i in 2:2) {
  
  rowGene = strsplit(line[i], "\t")[[1]] #contienetutti i valori del gene ordinati in modo decrescente
  
  nameGene = rowGene[1] #prendiamo il nome del gene
  
  #contienetutti i valori del gene ordinati in modo decrescente, eliminando il primo elemento
  #approssimazione 10 cifra dopo la virgola
  #le etichette sono all'interno di listsGene
  
  listsGene = sort.int(as.numeric(tail(rowGene,-1)), index.return = T,decreasing = T)
  
  #controllo se la riga considerata è tutta nulla
  if(any(listsGene[[1]] == 0) == FALSE) {
    lengthRow = lengths(listsGene[1]) #recupero la lunghezza della lista valuesGene
    
    valuesGene = unlist(listsGene[1])
    first_percentile = quantile(valuesGene, 0.33)
    second_percentile = quantile(valuesGene, 0.66)
    third_percentile = quantile(valuesGene, 0.99)
    
    div_array = split(valuesGene, ceiling(seq_along(valuesGene)/(ceiling(lengthRow/3))))
    
    #remove(valuesGene)
    
    up_array = as.numeric(unlist(div_array[1]))
    middle_array = as.numeric(unlist(div_array[2]))
    down_array = as.numeric(unlist(div_array[3]))
    
    #Se l'array è pari l'ultimo valore dell'ultimo array sarà la media, così la varianza non cambia
    
    if((lengthRow %% 2) == 1)
      down_array = c(down_array, mean(down_array))
    
    #remove(div_array)
    
    #Faccio la media dei tre gruppi
    
    mean_up = mean(up_array)
    mean_middle = mean(middle_array)
    mean_down = mean(down_array)
    
    #Fold change delle varie combinazioni (up vs mid, up vs down, etc..)
    if(dataIsLinear){
      
      fold_change_upvsmid = mean_middle/mean_up
      fold_change_upvsdown = mean_down/mean_up
      fold_change_midvsdown = mean_down/mean_middle
      
    }
    else{
      
      fold_change_upvsmid = 2^(abs(mean_middle-mean_up))
      fold_change_upvsdown = 2^(abs(mean_down-mean_up))
      fold_change_midvsdown =  2^(abs(mean_down-mean_middle))
      
    }
    
    #quali foldChange bisogna considerare che soddisfano il parametro dell'utente?
    
    if(fold_change_upvsmid > parametro && fold_change_upvsdown > parametro && fold_change_midvsdown > parametro){
      
      #differenza tra due array per evitare che il ttest fallisca
      if( (sd(middle_array-down_array) != 0) && (sd(up_array-middle_array) != 0) && (sd(up_array-down_array) != 0) ){
        
        ttest_up_mid1 = t.test(up_array, middle_array, var.equal=FALSE)
        ttest_up_down1 = t.test(up_array, down_array, var.equal=FALSE)
        ttest_mid_down1 = t.test(middle_array, down_array, var.equal=FALSE)
        #creare un data che contiene le etichette "cg.." che poi verranno composte alla fine del
        #json finale
        data <- list(
          list(
            ttest_up_mid = ttest_up_mid1$statistic[[1]],
            pvalue_up_mid =ttest_up_mid1$p.value,
            ttest_up_down = ttest_up_down1$statistic[[1]],
            pvalue_up_down = ttest_up_down1$p.value,
            ttest_mid_down = ttest_mid_down1$statistic[[1]],
            pvalue_mid_down = ttest_mid_down1$p.value,
            primo_percentile = first_percentile[[1]],
            secondo_percentile = second_percentile[[1]],
            terzo_percentile = third_percentile[[1]]))
        
        #vado alla ricerca delle etichette nel dataset GPL e 
        #succ. calcolerò i ttest trovati nel terzo dataset
        righeCheTiServono<- df[which(df$V22==nameGene), ]

        dataFrameFinale = NULL
        for(i in 1:length(righeCheTiServono$V1)){
          
          dataFrame<- df1[which(df1$sample==righeCheTiServono$V1[i]), ]
          #dataFrameFinale = rbind(dataFrameFinale,dataFrame[1, ])
        }
        
        dataFrameFinale = na.omit(dataFrameFinale)
        
  
        
        json1[[nameGene]] <- data
        
      } else {
        
        #creare un JSON con il nome del gene ed l'errore riportato
        paste("ERRORE nel calcolare il Ttest")
        errori=errori+1
        
        data <- list(
          list(
            errore = "ERRORE nel calcolare il Ttest"
          ))
        json1[[nameGene]] <- data
      }
    }
    else{
      data <- list(
        list(
          errore = "Il gene non rispetta il foldchange desiderato"
        ))
      json1[[nameGene]] <- data
    }
  }
  else{
    data <- list(
      list(
        errore = "Il gene ha tutti i valori settati a zero"
      ))
    json1[[nameGene]] <- data
  }
}


exportJson <- toJSON(json1)
write(exportJson, "test.json")
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
