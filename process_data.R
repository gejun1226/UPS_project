library(stringr)
library(data.table)
Sample_Info <- read.table("/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/Sample_Info.txt",header = T,stringsAsFactors = F,check.names = F)
RNA <- read.table("/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/RNAseq/RawData/ExprSet.txt",header = T,check.names = F)
Sample_Info <- Sample_Info[!duplicated(Sample_Info$ParticipantID),]
rownames(Sample_Info) <-Sample_Info$ParticipantID 
RNA_Info <- Sample_Info[colnames(RNA),]
RNA_Info2 <- RNA_Info[,7:11]
RNA <- colnames(RNA)
write.table(RNA_Info2,"/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/RNAseq/RawData/Group_Info.txt",sep = "\t",quote=F,row.names = T)
#RNA_sub <- RNA[,which(colnames(RNA) %in% Sample_Info$ParentSampleID)]
#a <- colnames(RNA_sub)
#name <- str_sub(a,1,9)
#colnames(RNA_sub) <- name
write.table(RNA_sub,"/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/RNAseq/RawData/ExprSet.txt",sep = "\t",quote=F)
Protein <- read.table("/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/Proteomics/RawData/ExprSet.txt",header = T,check.names = F)
Pro_Sample <- read.table("/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/Proteomics/RawData/Group_Info.txt",header = T)

Sample_Info <- Sample_Info[!duplicated(Sample_Info$AliquotID),]
rownames(Sample_Info) <- Sample_Info$AliquotID
#a <- Sample_Info[rownames(Pro_Sample),"ParticipantID"]
#rownames(Pro_Sample) <- a
write.table(Pro_Sample,"/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/Proteomics/RawData/Group_Info.txt",sep = "\t",quote=F,row.names = T)
RNA_Sample <- read.table("/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/RNAseq/RawData/Group_Info.txt",header = T)

rownames(Sample_Info) <- Sample_Info$ParticipantID
