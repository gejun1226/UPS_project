#========================KIRC
setwd("/Users/junge/Desktop/UPS/mass_data/KIRC/Sample_info")
KIRC <- read.table("S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.txt",sep="\t",stringsAsFactors = F,header = T)
KIRC_driver <- read.table("KIRC.txt",sep="\t",header = T,stringsAsFactors = F)
Tumor_KIRC <- KIRC[which(KIRC$Group == "Tumor"),]
Tumor_KIRC <- Tumor_KIRC[,c("ParticipantID","Parent.Sample.ID.s.","Aliquot.ID","TMT.plex","TMT.channel","Batch")]
colnames(Tumor_KIRC) <- c("ParticipantID","ParentSampleID","AliquotID","Plex","Channel","Batch")
KIRC_Sample_Info <- merge(Tumor_KIRC,KIRC_driver,by.x="ParticipantID",by.y="Tumor_Sample_Barcode")
write.table(KIRC_Sample_Info,"KIRC_Sample_Info_new.txt",sep = "\t",row.names = F,quote = F)
#=======================LUAD
setwd("/Users/junge/Desktop/UPS/mass_data/LUAD/Sample_info")
LUAD <- read.table("S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.txt",sep="\t",stringsAsFactors = F,header = T)
LUAD_driver <- read.table("LUAD.txt",sep="\t",header = T,stringsAsFactors = F)
Tumor_LUAD <- LUAD[which(LUAD$Type == "Tumor"),]
Tumor_LUAD <- Tumor_LUAD[,c("Participant.ID..case_id.","Broad.Sample.ID","Aliquot..Specimen.Label.","Experiment..TMT10plex.","Channel")]
Tumor_LUAD$Batch <- 0
colnames(Tumor_LUAD) <- c("ParticipantID","ParentSampleID","AliquotID","Plex","Channel","Batch")
LUAD_Sample_Info <- merge(Tumor_LUAD,LUAD_driver,by.x="ParticipantID",by.y="Tumor_Sample_Barcode")
write.table(LUAD_Sample_Info,"LUAD_Sample_Info_new.txt",sep = "\t",row.names = F,quote = F)
###################################################################################
#=======================UCEC==============================================
setwd("/Users/junge/Desktop/UPS/mass_data/UCEC_FBXW7/Sample_info/")
UCUC <- read.table("S043_CPTAC_UCEC_Discovery_Cohort_Study_Specimens_r1_Sept2018.txt",sep="\t",stringsAsFactors = F,header = T)
UCUC_driver <- read.table("UCEC.txt",sep="\t",stringsAsFactors = F,header = T)
Tumor_UCUC <- UCUC[which(UCUC$Group == "Tumor "),]
Tumor_UCUC <- Tumor_UCUC[,c("ParticipantID..Case_ID.","Parent.Sample.ID.s.","Aliquot.ID","TMT.plex","TMT.channel")]
Tumor_UCUC$Batch <- 0
colnames(Tumor_UCUC) <- c("ParticipantID","ParentSampleID","AliquotID","Plex","Channel","Batch")
#=====dot_index which will be 
dot_index <- grep(",",Tumor_UCUC$ParentSampleID)
#===========================
Tumor_UCUC_sub <- Tumor_UCUC[grep(",",Tumor_UCUC$ParentSampleID),]
ParentSampleID <- Tumor_UCUC_sub$ParentSampleID
ParentSampleID <- unlist(lapply(ParentSampleID, function(x) str_split(x,",")))
num_dot <- apply(Tumor_UCUC_sub,1,function(x) length(unlist(str_split(x["ParentSampleID"],","))))
ParticipantID <- rep(Tumor_UCUC_sub$ParticipantID,num_dot)
ParentSampleID <- ParentSampleID
AliquotID <- rep(Tumor_UCUC_sub$AliquotID,num_dot)
Plex <- rep(Tumor_UCUC_sub$Plex,num_dot)
Channel <- rep(Tumor_UCUC_sub$Channel,num_dot)
Batch <- rep(Tumor_UCUC_sub$Batch,num_dot)
Tumor_UCUC_sub <- data.frame(ParticipantID=ParticipantID,ParentSampleID=ParentSampleID,AliquotID=AliquotID,Plex=Plex,Channel=Channel,Batch=Batch)
#============================
Tumor_UCUC <- Tumor_UCUC[-dot_index,]
Tumor_UCUC_new<- rbind(Tumor_UCUC,Tumor_UCUC_sub)
Sample_Info <- merge(Tumor_UCUC_new,UCUC_driver,by.x="ParticipantID",by.y="Tumor_Sample_Barcode")
write.table(Sample_Info,"UCEC_Sample_Info_new.txt",sep = "\t",row.names = F,quote = F)
