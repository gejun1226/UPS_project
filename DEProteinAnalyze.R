#=====================================
packages <- c("stringr", "optparse", "limma","msmsEDA","msmsTests")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

lapply(packages, library, character.only = TRUE)
#=====================================================================
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="path to data", metavar="character"),
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="method of differential prteion calculation", metavar="character"),
  make_option(c("-e", "--model"), type="character", default=NULL, 
              help="Model of msmsTest, oviliable models are: msms.glm.pois, msms.glm.qlll, msms.edgeR.", 
              metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default='.', 
             help="path to output", metavar="character"),
  make_option(c("-n","--Na_Count_Cutoff"), type="character",default = NULL,help = "Set the NA count cutoff to select genes",metavar = "character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data_path <- opt$data_path
method <- opt$method
model <- opt$model
output_path <- opt$output_path
NaCutoff <- opt$Na_Count_Cutoff

#=========================
# DEProteinAnalyze function
#=========================
DEProteinAnalyze <- function(expr.path,SampleAnn.path,output.path,method="Limma",model=NULL,NaCutoff=0)
{
  exprSet <- read.table(expr.path,sep = "\t",header = T,stringsAsFactors = F,check.names = F)
  SampleAnn <- read.table(SampleAnn.path,sep = "\t",header = T,stringsAsFactors = F,check.names = F)
  genelist <- colnames(SampleAnn)
  
  na_count <- as.data.frame(apply(exprSet,1,function(x) sum(is.na(x))),stringsAsFactors = F)
  colnames(na_count) <- "count"
  numSam <- ncol(exprSet)
  #nacutoff <- as.numeric(round(numSam*NaCutoff))  # set the na cutoff at (sample)*0.3
  #na_count$count <- as.numeric(na_count$count)
  #selectGene <- rownames(na_count)[which(na_count$count <= nacutoff)]
  #exprSet <- exprSet[selectGene,]
  #exprSet_delNA <- exprSet[rownames(na_count)[which(na_count$count <= nacutoff)],] # expression which has values less than 30
  exprSet.c <- as.matrix(exprSet)
  exprSet.c[is.na(exprSet.c)] <- 0.000000001
  exprSet.c <- as.data.frame(exprSet.c) # mean_value
  sum <- as.data.frame(apply(exprSet.c, 1,function(x) sum(x)))  # calculate the sum of the data
  colnames(sum) <- "sum"
  sum_na <- merge(sum,na_count,by.x=0,by.y=0)
  sum_na$mean <- sum_na$sum/(numSam-sum_na$count)
  #contain_gene <- sum_na$Row.names[which(sum_na$count <= nacutoff)]
  #====================
  exprSet$mean <- sum_na$mean
  for ( i in seq(1:nrow(exprSet))){
    exprSet[i,which(is.na(exprSet[i,]))] <- exprSet[i,"mean"]  
  }
  #exprSet <- exprSet[contain_gene,]
  exprSet <- exprSet[,-ncol(exprSet)]
 #=====================
  for(i in genelist){
      # Make design matrix
      sub_sample <- SampleAnn[,i,drop=F]
      colnames(sub_sample) <- c("Condition")
      sub_sample$Condition[which((sub_sample$Condition)==0)] <- "Ctrl"  
      sub_sample$Condition[which((sub_sample$Condition)==1)] <- "Mutated"  
      sub_sample$Condition <-  factor(sub_sample$Condition,levels = c("Mutated","Ctrl"))
      exprSet <- exprSet[,rownames(sub_sample)]
      #=======================================
      if (method == "Limma"){
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
      } else if (method == "msmsTest"){
          if (!is.matrix(exprSet)){
              exprSet <- as.matrix(exprSet)
          }
          fd <- data.frame(otherfdata = rownames(exprSet))
          rownames(fd) <- rownames(exprSet)
          MSnSet_obj <- MSnSet(exprs=exprSet,fData=fd,pData=sub_sample)  # Note：pData=SampleAnn
          MSnSet_obj <- pp.msms.data(MSnSet_obj)  # pp.msms.data function used to deleted genes which all expression is 0.
          
          counts <- exprs(MSnSet_obj)
          null.f <- "y~1"
          alt.f <- "y~Condition"
          div <- apply(counts,2,sum)
          ### msmsTests method 
         if(model=="msms.glm.pois"){
              res <- msms.glm.pois(MSnSet_obj, alt.f, null.f, div=div) # msms.glm.pois；msms.glm.qlll；msms.edgeR
         }else if(model=="msms.glm.qlll"){
              res <- msms.glm.qlll(MSnSet_obj, alt.f, null.f, div=div)
         }else if(model=="msms.edgeR"){
              res <- msms.edgeR(MSnSet_obj, alt.f, null.f, div=div)
         }else{
              stop("Please check your model")
        }
        res$adj.p.value <-  p.adjust(res$p.value,method="BH")
        res <- res[,c("LogFC","p.value","adj.p.value","LR")]
        colnames(res) <- c("logFC","P.Value","adj.P.Val","LR")
    } else {
    stop("Method not available !!!")
    }
 write.table(res, paste(output.path,paste(i,".txt",sep = ""),sep = '/'),sep = "\t",quote = F)
    }
}
setwd("/Users/junge/Desktop/pepline_test/")
data_path <- "Data/SPOP/Proteomics"  # join(config["data"],"Proteomics")
output_path <- "Result/CPTAC/UCEC/Proteomics/DEP" # join(config["output"],"Proteomics/DEP")
method="Limma"
expr_path <- paste(data_path,"ExprSet.txt",sep = "/")
SampleAnn_path <- paste(data_path,"GroupInfo.txt",sep = "/")
dir.create(output_path,recursive = T)
DEProteinAnalyze(expr.path=expr_path,SampleAnn.path=SampleAnn_path,output.path = output_path,method=method,model=model)





