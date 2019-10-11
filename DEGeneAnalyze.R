#' Differential expression analysis
#' @docType methods
#' @name DEGeneAnalyze
#' @rdname DEGeneAnalyze
#' @param expr.path path to gene expression data
#' @param SampleAnn.path path to sample annotation data
#' The rownames of SampleAnn should match The colnames of exprSet, and the columns of SampleAnn should be the condition.
#' @param method Differential expression analysis method, e.g. limma, DESeq2.
#' @param output.path path to output
#' @return An ExprDataSet instance.

packages <- c("stringr", "optparse", "DESeq2","limma")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)

#===============================Parameter======================================
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="path to data", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default='.', 
              help="path to output", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='.', 
            help="", metavar="character"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data_path <- opt$data_path
output_path <- opt$output_path
method <- opt$method
#===============================DEGeneAnalyze======================================

DEGeneAnalyze <- function(expr.path,SampleAnn.path,output.path,method="DESeq2"){
  exprSet <- read.csv(expr.path,sep = "\t",stringsAsFactors = F,check.names = F)
  SampleAnn <- read.table(SampleAnn.path,sep = "\t",stringsAsFactors = F,check.names = F)
  genelist <- colnames(SampleAnn)
  #====================DESeq2================================
  if (method=="limma"){
    for(i in genelist){
      exprSet <- normalizeQuantiles(exprSet)
      sub_sample <- SampleAnn[,i,drop=F]
      colnames(sub_sample) <- c("Condition")
      sub_sample$Condition[which((sub_sample$Condition)==0)] <- "Ctrl"  
      sub_sample$Condition[which((sub_sample$Condition)==1)] <- "Mutated"  
      sub_sample$Condition <-  factor(sub_sample$Condition,levels = c("Mutated","Ctrl"))
      exprSet <- exprSet[,rownames(sub_sample)]
      
      design <- model.matrix(~0+sub_sample$Condition) # define matrix
      rownames(design) <- rownames(sub_sample)
      colnames(design) <- levels(sub_sample$Condition)
      # Make contrast matrix
      contrast.matrix<-makeContrasts(contrast="Mutated-Ctrl",levels = design) 
      # Limma model
      fit1 <- lmFit(exprSet, design=design) 
      fit2 <- contrasts.fit(fit1,contrasts = contrast.matrix)
      fit2 <- eBayes(fit2)
      res <- topTable(fit2, coef=1, n=Inf)
      res <- res[,c("logFC","P.Value","adj.P.Val","AveExpr","t","B")]
      write.table(res,paste(output.path,paste(i,".txt",sep=""),sep="/"),sep = "\t",quote = F)
    }
  } else if (method == "DESeq2"){
      for(i in genelist){
        sub_sample <- SampleAnn[,i,drop=F]
        colnames(sub_sample) <- c("Condition")
        sub_sample$Condition[which((sub_sample$Condition)==0)] <- "Ctrl"  
        sub_sample$Condition[which((sub_sample$Condition)==1)] <- "Mutated"  
        sub_sample$Condition <-  factor(sub_sample$Condition,levels = c("Ctrl","Mutated"))
        exprSet <- exprSet[,rownames(sub_sample)]
        colData <- data.frame(row.names = colnames(exprSet),condition=sub_sample$Condition)
        dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                      colData = colData,
                                      design = ~condition)
        dds_N <- DESeq(dds)
        #res <- results(dds_N)
        #res <- as.data.frame(res)
        resLFC <- lfcShrink(dds_N, coef="condition_Mutated_vs_Ctrl")
        resLFC <- as.data.frame(resLFC)
        resLFC <- resLFC[, c("log2FoldChange","pvalue","padj","baseMean","stat")]
        colnames(resLFC) <- c("logFC","P.Value","adj.P.Val","AveExpr","t")
        resLFC <- na.omit(resLFC)
        write.table(resLFC,paste(output.path,paste(i,".txt",sep=""),sep="/"),sep = "\t",quote = F)
      }
   }
}
setwd("/Users/junge/Desktop/pepline_test/")
data_path <- "Data/CPTAC/UCEC/RNAseq"  # join(config["input"],"RNAseq")
output_path <- "/Users/junge/Desktop/pepline_test/Result/CPTAC/UCEC/RNAseq/DEG" # join(config["output"],"RNAseq/DEP")
data_path <- "/Users/junge/Desktop/pepline_test/Data/SPOP/RNAseq"  # join(config["data"],"Proteomics")
output_path <- "/Users/junge/Desktop/pepline_test/Result/SPOP/RNAseq/DEG" # join(config["output"],"Proteomics/DEP")
method="Limma"
expr_path <- paste(data_path,"ExprSet.txt",sep = "/")
SampleAnn_path <- paste(data_path,"GroupInfo.txt",sep = "/")
dir.create(output_path,recursive = T)
DEGeneAnalyze(expr.path=expr_path,SampleAnn.path=SampleAnn_path,output.path = output_path,method=method)














