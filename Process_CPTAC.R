packages <- c("stringr", "optparse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, 
              help="path to data", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default='.', 
              help="path to output", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
data_path <- opt$data_path
output_path <- opt$output_path

####################
# functions
####################



Process_CPTAC <- function(expr.path, sample.info.path,output.path, selected="Unshared")
{     
      # Proteins Abundance Data
      expr <- read.delim(expr.path,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
      
      if(selected=="Unshared"){
      index <- grepl("Unshared.Log.Ratio",colnames(expr))
      }else{
      index <- !grepl("Unshared.Log.Ratio",colnames(expr))
      }
      
      sub_expr <- expr[,c(1,which(index))]
      rownames(sub_expr) <- sub_expr$Gene
      sub_expr$Gene <- NULL
      sub_expr <-sub_expr[-c(1:3),]
      # Sample information
      Sample_info <- read.table(sample.info.path,sep="\t",header = T,stringsAsFactors = F)
      
      # exprSet
      del_log_Ratio <- str_replace(colnames(sub_expr),"\\..*","")  #Delete the Log.Ration, change the samples name to Aliquot.ID
      Aliquot.ID <- Sample_info[,"AliquotID"]
      AliquotID_index <- grep("AliquotID",colnames(Sample_info))
      exprSet <- sub_expr[,which((del_log_Ratio) %in% Aliquot.ID)] 
      colnames(exprSet) <- str_replace(colnames(exprSet),"\\..*","") 
      exprSet_delna <- na.omit(exprSet)
      
      # MutationInfo
      #gene_index <- which(colnames(Sample_info) %in% mut_gene_list)
      gene_index <- seq(grep("Batch",colnames(Sample_info))+1,ncol(Sample_info))
      
      # Driver mutation sample
      Sample <- Sample_info[,c(AliquotID_index,gene_index)]
      Sample <- unique(Sample) # contain same sample names.
      rownames(Sample) <- Sample[,"AliquotID"] 
      Sample[,"AliquotID"] <- NULL
      
      # plex bacth sample
      Sample1 <- Sample_info[,-c(gene_index)]
      write.table(Sample1, paste(output.path,"SampleInfo2QC.txt",sep = "/"),sep = "\t",quote = F)
      write.table(exprSet, paste(output.path,"ExprSet.txt",sep = "/"),sep = "\t",quote = F)
     # write.table(exprSet_delna, paste(output.path,"ExprSet.txt",sep = "/"),sep = "\t",quote = F)
      write.table(Sample, paste(output.path,"Group_info.txt",sep = "/"),sep = "\t",quote = F)
}

### process CPTAC data
#data_path <- "/Users/junge/Desktop/MAUPS-master/Data/CPTAC" # config["data"]
#output_path <- "/Users/junge/Desktop/MAUPS-master/Result/CPTAC" # config["output"]

#data_path <- "/Users/junge/Desktop/MAUPS-master/Data/CPTAC" # config["data"]
#output_path <- "/Users/junge/Desktop/MAUPS-master/Result/CPTAC" # config["output"]

#i="UCEC"

for (i in list.files(data_path)){
  expr.path <- paste(data_path,i,'Protein_abundance.tsv',sep = '/')
  sample.info.path <- paste(data_path,i,'Sample_info.txt',sep = '/')
  output.path <- paste(output_path,i,"pre_process",sep = "/")
  dir.create(output.path,recursive = T)
  Process_CPTAC(expr.path=expr.path, sample.info.path=sample.info.path, selected="Unshared",output.path=output.path)
}

