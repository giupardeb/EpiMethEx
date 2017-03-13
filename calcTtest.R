start.time <- Sys.time()
library(RJSONIO)
#args <- commandArgs(TRUE)
conn <- file("genomicMatrix", open="r") #args dovr? contenere il path inviato da php
line <-readLines(conn)
close(conn)
errori=0
data=I
json1 <- list()
etichette = tail(strsplit(line[1],"\t")[[1]], -1)

for(i in 2:length(line)) {

rowGene = strsplit(line[i], "\t")[[1]] #contienetutti i valori del gene ordinati in modo decrescente

nameGene = rowGene[1] #prendiamo il nome del gene
#contienetutti i valori del gene ordinati in modo decrescente, eliminando il primo elemento
#approssimazione 10 cifra dopo la virgola

# valuesGene <- sort(as.numeric(tail(valuesGene,-1)), decreasing = T)

listsGene = sort.int(as.numeric(tail(rowGene,-1)), index.return=T,decreasing = T)
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
# #
# print("UP_ARRAY")
# print(up_array)
# # 
# print("////////////////////////////")
# print("MEAN_UP_ARRAY")
# print(mean(up_array))
# print("////////////////////////////")
# # 
# print("MIDDLE_ARRAY")
# print(middle_array)
# # 
# print("////////////////////////////")
# print("MEAN_MIDDLE_ARRAY")
# print(mean(middle_array))
# print("////////////////////////////")
# # 
# print("DOWN_ARRAY")
# print(down_array)
# # 
# print("////////////////////////////")
# print("MEAN_DOWN_ARRAY")
# print(mean(down_array))


#Se l'array è pari l'ultimo valore dell'ultimo array sarà la media, così la varianza non cambia

if((lengthRow %% 2) == 1)
  down_array = c(down_array, mean(down_array))

remove(div_array)

#differenza tra due array per evitare che il ttest fallisca

if( (sd(middle_array-down_array) != 0) && (sd(up_array-middle_array) != 0) && (sd(up_array-down_array) != 0) ){
  
  ttest_up_mid1 = t.test(up_array, middle_array, var.equal=FALSE)
  ttest_up_down1 = t.test(up_array, down_array, var.equal=FALSE)
  ttest_mid_down1 = t.test(middle_array, down_array, var.equal=FALSE)
  
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

# output = list(
#   t1=ttest_up_mid$statistic,
#   pvalue1=ttest_up_mid$p.value,
#   t2=ttest_up_down$statistic,
#   pvalue2=ttest_up_down$p.value,
#   t3=ttest_mid_down$statistic,
#   pvalue3=ttest_mid_down$p.value
# )


exportJson <- toJSON(json1)
write(exportJson, "test.json")
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
