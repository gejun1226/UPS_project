#=======================================================
packages <- c("optparse","ggplot2","ggrepel","stringr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)
option_list = list(
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="path to data", metavar="character"),
  make_option(c("-o", "--output_path"), type="character", default='.', help="path to output", metavar="character"),
  make_option(c("-p", "--pvalue_cutoff"), type="numeric", default=0.05, help="cutoff of pvalue", metavar="character"),
  make_option(c("-f", "--foldchange_cutoff"), type="numeric", default=1.5, help="cutoff of foldchange", metavar="character"),
  make_option(c("-l","--label"),type="character", default=NULL, help="labeled proteins or genes", metavar="character"),
  make_option(c("-t","--top"),type="character", default=10, help="top proteins or genes number", metavar="character"),
  make_option(c("-s","--select_gene_protein"),type="character", default="protein", help="select protein,gene or both", metavar="character")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#==========================================
data_path <- opt$data_path
output_path <- opt$output_path
pvalue_cut <- opt$pvalue_cutoff
FC_cut <- opt$foldchange_cutoff
label <- opt$label
top <- opt$top
select <- opt$select_gene_protein
#===========================================

VolcanoView <- function(expr.path, x = "logFC", y = "adj.P.Val", Label=NA, top = 0,
                        topnames = NULL, filename = NULL,
                        x_cutoff = 0.1, y_cutoff = 0.05, main = NULL,
                        xlab = "Log2 Fold Change", ylab = "-Log10(Adjust.P)", ...){
  mat <- read.table(expr.path,sep = "\t",header = T)
  gg = mat[,c(x,y)]
  gg$group="no"
  gg$group[gg[,x]>x_cutoff & gg[,y]<y_cutoff] = "up"
  gg$group[gg[,x]< -x_cutoff & gg[,y]<y_cutoff] = "down"
  gg[, y] = -log10(gg[, y])
  if(!(top==0 & is.null(topnames))){
    gg$Label = rownames(gg)
    if(!is.na(Label)) gg$Label = mat[, Label]
    gg = gg[order(gg[,y], decreasing = TRUE), ]
    idx1 = idx2 = c()
    if(top>0){
      idx1 = which(gg$group=="up")[1:top]
      idx2 = which(gg$group=="down")[1:top]
    }
    idx = unique(c(idx1, idx2, which(gg$Label %in% topnames)))
  }
  
  mycolour=c("no"="#969696",  "up"="#d7191c","down"="#2c7fb8","note"="#252525")
  #=========
  p = ggplot(gg, aes(x=gg[,x],y=gg[,y],colour=group,fill=group))
  p = p + geom_jitter(position = "jitter",show.legend = FALSE,size = 0.5)
  p = p + theme(text = element_text(colour="black",size = 14, family = "Helvetica"),
                plot.title = element_text(hjust = 0.5, size=16),
                axis.text = element_text(colour="black"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  p = p + geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted")
  p = p + geom_vline(xintercept = c(-x_cutoff, x_cutoff), linetype = "dotted")
  p = p + labs(x=xlab,y=ylab,title=main)
  p = p + annotate("text",color="black",x=x_cutoff, y=max(gg[,y]), hjust = 0, vjust = 1,
                   label=paste("Up: ",dim(gg[gg$group=="up",])[1],sep=""))
  p = p + annotate("text",color="black",x=(-x_cutoff), y=max(gg[,y]),hjust = 1, vjust = 1,
                   label=paste("Down: ",dim(gg[gg$group=="down",])[1],sep=""))
  if(!(top==0 & is.null(topnames)))
    p = p + ggrepel::geom_text_repel(aes(x=gg[idx,x],y=gg[idx,y], colour="note",label = Label), data=gg[idx,],
                                     fontface = 'bold', size = 2.5,
                                     box.padding = unit(0.5, "lines"), segment.color = 'black',
                                     point.padding = unit(0.5, "lines"), segment.size = 0.5)
  p = p + scale_color_manual(values=mycolour)
  p = p + scale_fill_manual(values=mycolour)
  p = p + theme(legend.position = "none")
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename,width = 5, height = 5)
  }
  return(p)
}

#=============================================
setwd("/Users/junge/Desktop/pepline_test/")
data_path  <- "Result/SPOP"  # join(config["output"])
output_path  <- "Result/SPOP/DEVolcanoView"  # join(config["output"],"DEVolcanoView")
select <- "Both"
FC_cut <- 1.5
pvalue_cut <- 0.05
top <- 10
label <- c("EP300","AR","ARRB1")
dir.create(output_path,recursive = T)
   if (select=="Both"){
    DEP_path  <- paste(data_path,"Proteomics","DEP",sep="/")
    for(i in list.files(DEP_path)){
    DEP.path <- paste(DEP_path,i,sep = "/")
    p_protein <- VolcanoView(expr.path=DEP.path, x = "logFC", y = "P.Value", Label=NA, top =top,
                           topnames = label, filename = paste(output_path,paste(str_replace_all(i,".txt",""),"Protein_VolcanoPlot",".pdf",sep=""),sep = "/"),
                           x_cutoff = log2(FC_cut), y_cutoff = pvalue_cut, main = NULL,
                           xlab = "Log2 Fold Change", ylab = "-Log10(P.Val)")
  }
    DEG_path  <- paste(data_path,"RNAseq","DEG",sep="/")
    for(i in list.files(DEG_path)){
    DEG.path <- paste(DEG_path,i,sep = "/") 
    p_gene <- VolcanoView(expr.path=DEG.path, x = "logFC", y = "P.Value", Label=NA, top = top,
                        topnames = label, filename = paste(output_path,paste(str_replace_all(i,".txt",""),"Gene_VolcanoPlot",".pdf",sep=""),sep = "/"),
                        x_cutoff = log2(FC_cut), y_cutoff = pvalue_cut, main = NULL,
                        xlab = "Log2 Fold Change", ylab = "-Log10(P.Val)",width = 5, height = 5)
     }
    }else if (select=="protein"){
    DEP_path  <- paste(data_path,"Proteomics","DEP",sep="/")
    for(i in list.files(DEP_path)){
      DEP.path <- paste(DEP_path,i,sep = "/")
      p_protein <- VolcanoView(expr.path=DEP.path, x = "logFC", y = "P.Value", Label=NA, top =top,
                               topnames = label, filename = paste(output_path,paste(str_replace_all(i,".txt",""),"Protein_VolcanoPlot",".pdf",sep=""),sep = "/"),
                               x_cutoff = log2(FC_cut), y_cutoff = pvalue_cut, main = NULL,
                               xlab = "Log2 Fold Change", ylab = "-Log10(P.Val)")
    }
    } else if (select=="gene"){
      DEG_path  <- paste(data_path,"RNAseq","DEG",sep="/")
      for(i in list.files(DEG_path)){
        DEG.path <- paste(DEG_path,i,sep = "/") 
        p_gene <- VolcanoView(expr.path=DEG.path, x = "logFC", y = "P.Value", Label=NA, top = top,
                              topnames = label, filename = paste(output_path,paste(str_replace_all(i,".txt",""),"Gene_VolcanoPlot",".pdf",sep=""),sep = "/"),
                              x_cutoff = log2(FC_cut), y_cutoff = pvalue_cut, main = NULL,
                              xlab = "Log2 Fold Change", ylab = "-Log10(P.Val)",width = 5, height = 5)
      }
  }

