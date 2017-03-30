start.time <- Sys.time()
library(data.table)
df = fread("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/probe/GPL13534-11288.txt", sep="\t",header=F,skip=37)
righeCheTiServono<- df[which(df$V22=="MMP9"), ]

df1 = fread("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/Methylation450k/genomicMatrix", sep="\t")
#dataFrameFinale = data.frame(subset(df1, FALSE))
dataFrameFinale = NULL
for(i in 1:length(righeCheTiServono$V1)){
  
dataFrame<- df1[which(df1$sample==righeCheTiServono$V1[i]), ]
dataFrameFinale = rbind(dataFrameFinale,dataFrame[1, ])
}

dataFrameFinale = na.omit(dataFrameFinale)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
