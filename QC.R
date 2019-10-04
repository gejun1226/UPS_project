packages <- c("ggplot2","stringr", "optparse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

lapply(packages, library, character.only = TRUE)
#=====================================================================
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="path to data", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default='.', 
              help="path to output", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_path <- opt$data_path
output_path <- opt$output_path

QC <- function(expr.path,sample.info.path,output.path){
    exprSet <- read.table(expr.path,sep = "\t",header = T,stringsAsFactors = F,check.names = F)
    #  sample.info
    na_count <- as.data.frame(apply(exprSet,1,function(x) sum(is.na(x))),stringsAsFactors = F)
    colnames(na_count) <- "count"
    #ggplot(na_count,aes(x=count))+geom_density()
    count_percent <- as.data.frame(table(na_count),stringsAsFactors = F)
    count_percent$Ratio <- round(count_percent$Freq/nrow(na_count),4)
    count_percent$na_count <- as.integer(unlist(count_percent$na_count))
    Sample <-  ncol(exprSet)
    Gene <- nrow(exprSet)
    Num0 <- count_percent[which(count_percent$na_count==0),"Freq"]
    count_percent <- count_percent[-1,]
    count_percent$sum <- count_percent$Ratio[1]
    for(i in seq(2,nrow(count_percent))){
      count_percent$sum[i] <- count_percent$sum[i-1]+count_percent$Ratio[i]
    }
    p_NA <- ggplot(count_percent, aes(x=na_count,y=Ratio)) + geom_histogram(stat = "identity",bins=1,fill="#4292c6")+theme_bw()
     
    p_NA <- geom_density(data=na_count,mapping=aes(x=count))+p_NA
    
    p_NA <- p_NA + geom_text(aes(label=as.character(Freq)),vjust=-0.5,size = 2.5)
    geom_density(mapping = NULL, data = NULL, stat = "density",
                 position = "identity", ..., na.rm = FALSE, show.legend = NA,
                 inherit.aes = TRUE)
    #p_NA <- p_NA + geom_line(count_percent,aes(x=as.factor(na_count),y=sum), stat="identity", colour="skyblue", size = 1)
    
    #p <- p + geom_text(aes(label=as.character(round(sum,2)),vjust=-0.8))
    p_NA <- p_NA + theme(axis.text.x = element_text(size =8,angle=60,face = "bold"))
    p_NA <- p_NA + theme(axis.text.y = element_text(size = 8,face="bold"))
    p_NA <- p_NA + xlab("")+ylab("")
    p_NA <- p_NA + theme(plot.title = element_text(face = "bold",size=15))
    p_NA <- p_NA + theme(
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank()
        
    )
    p_NA <- p_NA + annotate("text",x=max(count_percent$na_count),y=max(count_percent$Ratio),label=paste("Samples: ",Sample,sep = ""),size = 4,hjust = 3, vjust = 1)
    p_NA <- p_NA + annotate("text",x=max(count_percent$na_count),y=max(count_percent$Ratio),label=paste("Total Genes: ",Gene,sep = ""),size = 4,hjust = 2, vjust = 3)
    p_NA <- p_NA + annotate("text",x=max(count_percent$na_count),y=max(count_percent$Ratio),label=paste("No NAvalue Genes: ",Num0,sep = ""),size = 4,hjust = 1.5, vjust = 5)
    ggsave(paste(output.path,"NADistribution.pdf.pdf",sep = "/"),p_NA,height = 5,width=10)

   #=================================================== 
    SampleInfo <- read.table(sample.info.path,sep = "\t",header = T,stringsAsFactors = F,check.names=F)
    identified_protein <- as.data.frame(nrow(exprSet)-apply(exprSet,2,function(x) sum(is.na(x))))
    colnames(identified_protein) <- c("Identified.Count")
    info <- colnames(SampleInfo)
  
    #===================各个样本的count信息
    identified_protein$Sample<- as.factor(rownames(identified_protein))
    p <- ggplot(data=identified_protein,aes(x=Sample,y=Identified.Count,colour=Sample))+geom_point()+theme_bw()
    p <- p + theme(legend.position = "none") 
    
    p <- p + ylim(min(identified_protein$Identified.Count)-2000,max(identified_protein$Identified.Count)) 
    p <- p + labs(x="",y="",title = "Identified Proteins")
    p <- p + theme(plot.title = element_text(hjust = 0.5,face = "bold",size=15)) # 调整标题
    p <- p + theme(axis.text.x = element_text(size = 8, angle = 90,face = "bold")) # 调整x轴
    
    p <- p + theme(axis.text.y = element_text(size = 10,face = "bold"))  # 调整y轴
    p <- p + theme(panel.border = element_blank()) # 去掉外边框线
    p <- p + theme(axis.line = element_line(size=0.3, colour = "black")) # 加上坐标轴
   
    ggsave(paste(output.path,"IdentifiedPro.pdf",sep = "/"),p,height = 5,width=10)
    #===================
    exprSet_delna <- na.omit(exprSet)
    exprSet_delna_pca <- prcomp(t(exprSet_delna))
    pca_out <- as.data.frame(exprSet_delna_pca$x)
    pca_out$Channel <- as.factor(SampleInfo[rownames(SampleInfo),"Channel"])
    
    pca_out$Plex <- as.factor(SampleInfo[rownames(SampleInfo),"Plex"])
    pca_out$Batch <- as.factor(SampleInfo[rownames(SampleInfo),"Batch"])
    percentage <- round(exprSet_delna_pca$sdev / sum(exprSet_delna_pca$sdev) * 100, 2) #取百分比
    percentage <- paste(colnames(pca_out), "(", as.character(percentage), "%", ")", sep="")
    
    # Channel Batch
    if("Channel" %in% info){
        p_Channel <- ggplot(pca_out,aes(x=PC1,y=PC2,color=Channel))+theme_bw()
        p_Channel <- p_Channel+geom_point()+ xlab(percentage[1]) + ylab(percentage[2])
        ggsave(paste(output.path,"PCA_Channel.pdf",sep = "/"),p_Channel,width = 6,height = 4)
     }
    # Plex Batch
    if("Plex" %in% info){
        p_Plex <- ggplot(pca_out,aes(x=PC1,y=PC2,color=Plex))+theme_bw()
        p_Plex <- p_Plex+geom_point()+ xlab(percentage[1]) + ylab(percentage[2])
        ggsave(paste(output.path,"PCA_Plex.pdf",sep = "/"),p_Plex,width = 6,height = 5)
    }
    # Batch
    if("Batch" %in% info){
        p_Batch <- ggplot(pca_out,aes(x=PC1,y=PC2,color=Batch))+theme_bw()
        p_Batch <- p_Batch+geom_point()+ xlab(percentage[1]) + ylab(percentage[2])
        ggsave(paste(output.path,"PCA_Batch.pdf",sep = "/"),p_Batch,width = 6,height = 5)
    }
}
data_path <- "/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/Proteomics"  # join(config["input"],"Proteomics")
#output_path <- "/Users/junge/Desktop/pepline_test/Data/CPTAC/UCEC/Proteomics/QC" # join(config["output"],"Proteomics/QC")
expr.path <- paste(data_path,"ExprSet.txt",sep = "/")
sample.info.path <- paste(data_path,"Batch_Info.txt",sep = "/")
dir.create(output_path,recursive = T)
QC(expr.path=expr.path,sample.info.path=sample.info.path,output.path=output_path)


