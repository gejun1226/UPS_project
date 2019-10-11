packages <- c("stringr", "optparse","VennDiagram")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)
#===============================Parameter======================================
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="path to data", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default='.', help="path to output", metavar="character"),
  make_option(c("-p", "--pvalue_cutoff"), type="numeric", default=0.05, help="cutoff of pvalue", metavar="character"),
  make_option(c("-f", "--foldchange_cutoff"), type="numeric", default=1.5, help="cutoff of foldchange", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
data_path <- opt$data_path
output_path <- opt$output_path
pvalue_cut <- opt$pvalue_cutoff
FC_cut <- opt$foldchange_cutoff

OverlapDEProGene <- function(DEG.path,DEP.path,output.path,pvalue_cut=0.05,FC_cut=1.5){
   DEG <- read.table(DEG.path,sep = "\t")
   DEP <- read.table(DEP.path,sep = "\t")
   overlap_gene <- intersect(rownames(DEG),rownames(DEP))
   com_gene <- length(overlap_gene)

   DEG <- DEG[overlap_gene,]
   DEP <- DEP[overlap_gene,]
   
   up_gene <- DEG[DEG$P.Value < pvalue_cut & DEG$logFC > log2(FC_cut),]
   up_protein <- DEP[DEP$P.Value < pvalue_cut & DEP$logFC > log2(FC_cut),]
   down_gene <- DEG[DEG$P.Value < pvalue_cut & DEG$logFC < -log2(FC_cut),]
   down_protein <- DEP[DEP$P.Value < pvalue_cut & DEP$logFC < -log2(FC_cut),]
   venn.diagram(x=list(up_gene=rownames(up_gene),up_protein=rownames(up_protein),down_gene=rownames(down_gene),down_protein=rownames(down_protein)),
                col = c("#238b45", "#fdd49e", "#8c96c6", "#08519c"),fill=c("#238b45", "#fdd49e", "#8c96c6", "#08519c"),
                filename = paste(output.path,paste(str_replace_all(i,".txt",""),"_Venn",".pdf",sep=""),sep = "/"))
}
setwd("/Users/junge/Desktop/pepline_test/")
data_path  <- "Result/SPOP"  # join(config["output"])
output_path  <- "Result/SPOP/CompaProRNA"  # join(config["output"],"CompaProRNA")
DEG_path  <- paste(data_path,"RNAseq","DEG",sep="/")
DEP_path  <- paste(data_path,"Proteomics","DEP",sep="/")
dir.create(output_path,recursive = T)
for(i in list.files(DEG_path)){
  DEG.path <- paste(DEG_path,i,sep = "/") 
  DEP.path <- paste(DEP_path,i,sep = "/")
  OverlapDEProGene(DEG.path=DEG.path,DEP.path=DEP.path,output.path=output_path,pvalue_cut=pvalue_cut,FC_cut=FC_cut)
}





