#' GetCptacRNAseqData: Together the RNAseq data(Count,FPKM) from NCI Genomic Data Commons (GDC)
#' @description : UCEC,Kidney and Lung cancer RNAseq data are stored in the NCI Genomic Data Commons (GDC).
#'               The data which we get from the GDC is separate.We can use GetCptacRNAseqData function
#'               to merge the data into one set and transform the "ensembl_gene_id_version" to "hgnc_symbol".
#' @param folder_name: The name of your folder，where stores the count(or FPKM) data downloaded from GDC.
#' @param annofile_name: The name of the comment file which is downloaded from GDC.
#' @param path: The path to the all files.
#' @param type: Data type (Count or FPKM).
#' @return Merged data, which can be used to do the downstream analysis.
GetCptacRNAseqData <- function(folder_name="gdc_download", annofile_name = "gdc_sample_sheet.tsv",path="...",type="FPKM"){
  
      setwd(path)
      folder_path <- paste(path,folder_name,sep = "/")
      dirs <- list.files(folder_path)
      dirs <- dirs[!grepl(".txt",dirs)]
      # file_id
      anno_file_path <- paste(path,annofile_name,sep = "/")
      anno_file <- read.table(anno_file_path,sep = "\t",header = T,stringsAsFactors = F)
      rownames(anno_file) <- anno_file[,1]
      anno_file <- anno_file[,-1]
      anno_file <- anno_file[dirs,]
      files <- paste(folder_name,dirs,anno_file$File.Name,sep = "/")
      
      
      if (type=="FPKM"){
          df <-  lapply(files,function(x) read.csv(file=x, sep = "\t",stringsAsFactors = F,header=F)) %>% bind_cols() 
          gene_expr_which <- c(1,seq(2,ncol(df),2))
          gene_expr <- df[,gene_expr_which]
      }else {
          df <-  lapply(files,function(x) read.csv(file=x, sep = "\t",stringsAsFactors = F,header=T)) %>% bind_cols()
          df <- df[-c(1:4),]
          gene_expr_which <- c(1,which(grepl("unstranded",colnames(df))))
          gene_expr <- df[,gene_expr_which]
      }  
      colnames(gene_expr) <- c("entrze",anno_file$Sample.ID)
      colnames(gene_expr) <- str_replace_all(colnames(gene_expr),",.*","")
      gene_expr$entrze <- str_replace_all(gene_expr$entrze,"\\..*","")
      
      
      # Use biomart to transform gene ID.
      requireNamespace("biomaRt", quietly=TRUE) || stop("need biomaRt package")
          ensembl <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
          filters <-  listFilters(ensembl)      #ensembl_gene_id_version   # hgnc_symbol
          ense_sym <-getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = gene_expr$entrze,
                   mart = ensembl)
          ense_sym <- ense_sym[!ense_sym$hgnc_symbol=="",]
          gene_expr_symbol <- merge(gene_expr,ense_sym,by.x=1,by.y=1)
          gene_expr_unique <- gene_expr_symbol[!duplicated(gene_expr_symbol$hgnc_symbol),]
          rownames(gene_expr_unique) <- gene_expr_unique[,"hgnc_symbol"]
          gene_expr_unique <- gene_expr_unique[,-as.numeric(ncol(gene_expr_unique))]
          gene_expr_unique <- gene_expr_unique[,colnames(gene_expr_unique)!="entrze"]
      #write.table(gene_expr_unique,"CPTAC_UCEC_RNAseq_FPKM.txt",sep = "\t",quote = F)
      return(gene_expr_unique)
      }
