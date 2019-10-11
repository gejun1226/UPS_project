#' Correlation between RNAseq and Protein
#' @docType methods
#' @name CorrRNAPro
#' @rdname CorrRNAPro
#' @param Pro.path Path to protein abundance data
#' @param RNA.path Path to RNAseq data
#' @param output.path 
#' @return Corrlation plot.
#' @importFrom import corr.test from psych
#=======================================================
packages <- c("optparse","ggplot2","psych") # import corr.test from psych
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="path to data", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default=NULL, 
              help="path to output", metavar="character"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data_path <- opt$data_path
output_path <- opt$output_path
#========================================================
      DataCorrelation <- function(Pro.path,RNA.path,output.path){
          Protein <- read.table(Pro.path,sep = "\t",stringsAsFactors = F,check.names = F) #contain NA  note:should use the data which fill na
          #============================================
          na_count <- as.data.frame(apply(Protein,1,function(x) sum(is.na(x))),stringsAsFactors = F)
          colnames(na_count) <- "count"
          numSam <- ncol(Protein)
          Protein.c <- as.matrix(Protein)
          Protein.c[is.na(Protein.c)] <- 0.000000001
          Protein.c <- as.data.frame(Protein.c) # mean_value
          sum <- as.data.frame(apply(Protein.c, 1,function(x) sum(x)))  # calculate the sum of the data
          colnames(sum) <- "sum"
          sum_na <- merge(sum,na_count,by.x=0,by.y=0)
          sum_na$mean <- sum_na$sum/(numSam-sum_na$count)
          #contain_gene <- sum_na$Row.names[which(sum_na$count <= nacutoff)]
          #====================
          Protein$mean <- sum_na$mean
          for ( i in seq(1:nrow(Protein))){
            Protein[i,which(is.na(Protein[i,]))] <- Protein[i,"mean"]  
          }
          #Protein <- Protein[contain_gene,]
          Protein <- Protein[,-ncol(Protein)]
          #====================
          #Protein <- na.omit(Protein)
          RNA <- read.table(RNA.path,sep = "\t",stringsAsFactors = F,check.names = F)
          overlap_gene <- intersect(rownames(Protein),rownames(RNA))
          com_gene <- length(overlap_gene)
          overlap_sample <- intersect(colnames(Protein),colnames(RNA))
          com_sample <- length(overlap_sample)
          Protein_sub <- Protein[overlap_gene,overlap_sample]
          RNA_sub <- RNA[overlap_gene,overlap_sample]
          #=========================Sample correlation
          p <- c()
          r <- c()
          for(i in seq(1,com_sample)){
              Protein_test <- Protein_sub[,i,drop=F]
              RNA_test <- RNA_sub[,i,drop=F]
              cor <- corr.test(Protein_test,RNA_test,method="spearman")
              p <- rbind(p,cor$p)
              r <- rbind(r,cor$r)
          }
          res <- as.data.frame(cbind(p,r))
          colnames(res) <- c("p.value","corr")
          p <- ggplot(res, aes(x=corr)) + geom_histogram(colour="#000000", fill="#0868ac")+theme_bw()  # Overlay with transparent density plot
          p <- p + theme(axis.line = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.border = element_blank(), panel.background = element_blank())
          p <- p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                        plot.title = element_text(hjust = 0.5, size=16),
                        axis.text = element_text(colour="black"))
          p <- p + theme(axis.text.x = element_text(size = 10,face = "bold"))
          p <- p + theme(axis.text.y = element_text(size = 10,face="bold"))
          p <- p + xlab("Correlation")
          p <- p + annotate("text",x=max(res$corr),y=com_sample/10,label=paste("Common_Samples: ",com_sample,sep = ""),size = 4,hjust = 3, vjust = 1)
          #========== significant ================
          # res$si <- ""
          # res$si[which(res$p.value <0.01)] <- "***"
          # res$si[which(res$p.value >=0.01 & res$p.value <0.05)] <- "**"
          # res$si[which(res$p.value >= 0.05 & res$p.value < 0.1)] <- "*"
          # p <- ggplot(res, aes(x=as.factor(rownames(res)),y=corr)) + geom_histogram(stat = "identity",bins=1,fill="#4292c6")+theme_bw()
          # p <- p + geom_text(aes(label=si),vjust=-0.4,size = 2.5)
          # p <- p + theme(axis.text.x = element_text(size =8,angle=90,face = "bold"))
          # p <- p + theme(axis.text.y = element_text(size = 8,face="bold"))
          # p <- p + xlab("")+ylab("Correlation ")
          # p <- p + theme(axis.title.y = element_text(face = "bold",size=15))
          # p <- p + theme(
          #   panel.grid.major = element_blank(),
          #   panel.grid.minor = element_blank(),
          #   panel.background = element_blank()
          # )
          # p <- p + annotate("text",x=nrow(res),y=max(res$corr),label=paste("Common Samples: ",com_sample,sep = ""),size = 3,hjust = 1, vjust = 1)
          # p <- p + annotate("text",x=nrow(res),y=max(res$corr),label=paste("Common Genes: ",com_gene,sep = ""),size = 3,hjust = 1, vjust = 3)
          ggsave(paste(output.path,"Sample_Correlation.pdf",sep = "/"),p,height = 4,width=5)
          #================== Gene Correlation
          p_gene <- c()
          r_gene <- c()
          for(i in seq(1,com_gene)){
            Protein_test <- t(Protein_sub[i,,drop=F])
            RNA_test <- t(RNA_sub[i,,drop=F])
            cor <- corr.test(Protein_test,RNA_test,method="spearman")
            p_gene <- rbind(p_gene,cor$p)
            r_gene <- rbind(r_gene,cor$r)
          }
          res_gene <- as.data.frame(cbind(p_gene,r_gene))
          res_gene <- na.omit(res_gene)
          colnames(res_gene) <- c("p.value","corr")
          pg <- ggplot(res_gene, aes(x=corr)) + geom_histogram(colour="#000000", fill="#0868ac")+theme_bw()  # Overlay with transparent density plot
          pg <- pg + theme(axis.line = element_line(size=0.5, colour = "black"),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.border = element_blank(), panel.background = element_blank())
          pg <- pg + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                         plot.title = element_text(hjust = 0.5, size=16),
                         axis.text = element_text(colour="black"))
          pg <- pg + theme(axis.text.x = element_text(size = 10,face = "bold"))
          pg <- pg + theme(axis.text.y = element_text(size = 10,face="bold"))
          pg <- pg + xlab("Correlation")
          pg <- pg + geom_density(alpha=.2)
          pg <- pg + annotate("text",x=max(res_gene$corr),y=com_gene/10,label=paste("common_gene: ",com_gene,sep = ""),size = 4,hjust = 3, vjust = 1)
          ggsave(paste(output.path,"Gene_Correlation.pdf",sep = "/"),pg,height = 4,width=5)
      }

#================================================================================
setwd("/Users/junge/Desktop/pepline_test/")
data_path <- "Data/SPOP"  # join(config["input"])    
output_path <- "Result/SPOP/DataCorrelation"  # join(config["output"],"DataCorrelation")
Pro_path <- paste(data_path,"Proteomics","ExprSet.txt",sep = "/")   
RNA_path <- paste(data_path,"RNAseq","ExprSet.txt",sep = "/") 
dir.create(output_path,recursive = T)
DataCorrelation(Pro.path=Pro_path,RNA.path=RNA_path,output.path = output_path)

