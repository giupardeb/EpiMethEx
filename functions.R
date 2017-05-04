metilazione <- function() {
  i <<- 1
  data <<- list()
  
  switch(c('a','b','c','d','e','f','g')[fc_Test],
         "a"={
           #upvsmid
           data <- ttester("UPvsMID",rowGeneUP,rowGeneMID)
         },
         "b"={
           #upvsdown
           data <- ttester("UPvsDOWN",rowGeneUP,rowGeneDOWN)
         },
         "d"={
           #midvsdown
           data <- ttester("MIDvsDOWN",rowGeneMID,rowGeneDOWN)
         },
         "c"={
           #upvsmid
           data <- ttester("UPvsMID",rowGeneUP,rowGeneMID)
           #upvsdown
           data <- ttester("UPvsDOWN",rowGeneUP,rowGeneDOWN)
         },
         "f"={
           #upvsdown
           data <- ttester("UPvsDOWN",rowGeneUP,rowGeneDOWN)
           #midvsdown
           data <- ttester("MIDvsDOWN",rowGeneMID,rowGeneDOWN)
         },
         "e"={
           #upvsmid
           data <- ttester("UPvsMID",rowGeneUP,rowGeneMID)
           #midvsdown
           data <- ttester("MIDvsDOWN",rowGeneMID,rowGeneDOWN)
         },
         {
           #prendi tutti e tre gli array
           data <- ttester("UPvsMID",rowGeneUP,rowGeneMID)
           data <- ttester("UPvsDOWN",rowGeneUP,rowGeneDOWN)
           data <- ttester("MIDvsDOWN",rowGeneMID,rowGeneDOWN)
         })
  
  return(data)
}

ttester <- function(name,array1,array2){
  if((sd(array1-array2) != 0)){
    
    ttester = t.test(array1, array2, var.equal=FALSE)
    
    data[[i]] <- paste(paste("ttest",name,sep = "_"),ttester$statistic[[1]],sep = "=")
    i<-i+1
    data[[i]] <- paste(paste("pvalue",name,sep = "_"),ttester$p.value,sep = "=")
      
  }
  else{
    #creare un JSON con il nome del gene ed l'errore riportato
    errori = errori + 1
    
    data[[i]] <- paste("errore = ERRORE nel calcolare il Ttest")
    i<-i+1
  }
  return(data)
}

splitdf <- function(df, n) {
  indx <- matrix(seq_len(ncol(df)), ncol = n)
  lapply(seq_len(n), function(x) df[, indx[, x]])
}