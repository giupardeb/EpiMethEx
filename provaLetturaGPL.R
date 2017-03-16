df = read.table(pipe("cut -f1,22 /home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/probe/GPL13534-11288.txt"),sep="\t",header=F,skip=37)
righeCheTiServono<- df[which(df$V2=="MMP9"), ]

df1 <- read.table(pipe("cut -f1 /home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/Methylation450k/genomicMatrix"),sep="\t")

#DataFrameFinale <- df1[which(df1, match(df1$V1,righeCheTiServono$V1,nomatch = FALSE)), ]
DataFrameFinale <- which(is.element(df1$V1,righeCheTiServono$V1))
for(i in DataFrameFinale) {
  prova <- read.table("/home/giupardeb/Scrivania/TesiRepoGithub/utilitiesTesi/dataset/Methylation450k/genomicMatrix",nrows=1, skip=i-1, sep="\t")
  
}
