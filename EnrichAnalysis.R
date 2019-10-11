#' Enrichment analysis

packages <- c("stringr", "optparse","MAGeCKFlute")   # import EnrichAnalyzer & EnrichedView from MAGeCKFlute
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)
#================================================
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="path to data", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default='.', help="path to output", metavar="character"),
  make_option(c("-s","--select_gene_protein"),type="character", default="protein", help="select protein,gene or both to do the enrichment analysis", metavar="character"),
  make_option(c("-m","--method_ORT_GSEA_HGT"),type="character",default = "GSEA",help = "method to do the enrichment analysis",metavar = "character"),
  make_option(c("-t","--top pathway"),type="numeric",default = 5,help = "top pathway to show",metavar = "character"),
  make_option(c("-b","--bottom pathway"),type="numeric",default = 5,help = "bottom pathway to show",metavar = "character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data_path <- opt$data_path
output_path <- opt$output_path
select <- opt$select_gene_protein
method <- opt$method_ORT_GSEA_HGT
top <- opt$top
bottom <- opt$bottom
#=================================================
Enrich <- function(DE.path,output.path,type="protein",method="GSEA",top=5,bottom=5){
  DE <- read.table(DE.path,sep = "\t")
  index <- grepl("log",colnames(DE))
  geneList <- DE[,index]
  names(geneList) <- rownames(DE)
  enrich <- EnrichAnalyzer(geneList = geneList, keytype = "Symbol", 
                          method = method, type = "GOBP+KEGG+CORUM")
  EnrichedView(enrich, rank_by = "p.adjust", top = top, bottom = bottom,
               x = "LogP", filename = paste(output.path,paste(str_replace(i,".txt",""),method,type,"enrich.pdf",sep = "_"),sep = "/"))
}
setwd("/Users/junge/Desktop/pepline_test/")
data_path  <- "Result/SPOP"  # join(config["output"])
output_path  <- "Result/SPOP/EnrichResult"  # join(config["output"],"Pathway")
#i = "CTNNB1.txt"
method="ORT"
type="Both"
top=20
bottom=0
dir.create(output_path,recursive = T)
if (select=="Both"){
  DEG_path  <- paste(data_path,"RNAseq","DEG",sep="/")
  DEP_path  <- paste(data_path,"Proteomics","DEP",sep="/")
  for(i in list.files(DEP_path)){
    DEP.path <- paste(DEP_path,i,sep = "/")
    Enrich(DE.path=DEP.path,output.path=output_path,type="protein",method = method,top=top,bottom = bottom)
    DEG.path <- paste(DEG_path,i,sep = "/")
    Enrich(DE.path=DEG.path,output.path=output_path,type="gene",method = method,top=top,bottom = bottom)
  }
  }else if (select=="protein"){
    DEP_path  <- paste(data_path,"Proteomics","DEP",sep="/")
    for(i in list.files(DEP_path)){
      DEP.path <- paste(DEP_path,i,sep = "/")
      Enrich(DE.path=DEP.path,output.path=output_path,type="protein",method = method,top=top,bottom = bottom)
    }
  } else if (select=="gene"){
    DEG_path  <- paste(data_path,"RNAseq","DEG",sep="/")
    for(i in list.files(DEG_path)){
      DEG.path <- paste(DEG_path,i,sep = "/")
      Enrich(DE.path=DEG.path,output.path=output_path,type="gene",method = method,top=top,bottom = bottom)
    }
  }






