#=======================================================
packages <- c("optparse","ggplot2","biomaRt") 
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
#==================================================================================
Rabit_input <- function(expr.path){
  ensembl <- useMart("ensembl")
  ensembl <-  useMart("ensembl",dataset = "hsapiens_gene_ensembl")
  diff_gene <- read.table(expr.path,sep = "\t",stringsAsFactors = F)
  TransId <- getBM(attributes=c('hgnc_symbol', 'entrezgene_id'), 
        filters = 'hgnc_symbol', 
        values = rownames(diff_gene), 
        mart = ensembl)
  TransId <- na.omit(TransId)
  diff_gene <- merge(diff_gene,TransId,by.x = 0,by.y=1)
  diff_gene <- diff_gene[!duplicated(diff_gene$entrezgene_id),]
  rownames(diff_gene) <- diff_gene[,"entrezgene_id"]
  diff_gene <- diff_gene[,!names(diff_gene) %in% c("Row.names", "entrezgene_id")]
  index <-  grepl("log",colnames(diff_gene))
  rabit_input <- diff_gene[,index,drop=F]
  write.table(rabit_input,paste(output_path,i,sep="/"),sep = "\t",quote = F,row.names = T) #diff_gene output.path新建一个output:Rabit
}
setwd("/Users/junge/Desktop/pepline_test/")
data_path  <- "Result/SPOP"  # join(config["output"])
output_path  <- "Result/SPOP/Rabit"  # join(config["output"],"Rabit")
DEG_path  <- paste(data_path,"RNAseq","DEG",sep="/")
dir.create(output_path,recursive = T)
for(i in list.files(DEG_path)){
  expr_path <- paste(DEG_path,i,sep = "/") 
  Rabit_input(expr.path=expr_path)
}

