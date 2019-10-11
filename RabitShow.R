#=======================================================
packages <- c("optparse","ggplot2") 
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
#=======================================================
# Rabit to show

RabitView <- function(rabitres.path,output.path,n=20){
    rabit <- read.table(rabitres.path,sep = "\t",header = T,stringsAsFactors = F)
    colnames(rabit) <- c("TF","Cp","Estimate","SE","t.value","p.value")
    rabit <- rabit[-c(1:3),]
    rabit$TF <- str_extract(rabit$TF,"\\..*@")
    rabit$TF  <- str_replace(rabit$TF,".","")
    rabit$TF  <- str_replace(rabit$TF,"@","")
    rabit <- rabit[order(rabit$p.value,decreasing = FALSE),]
    rabit <- rabit[!duplicated(rabit$TF),]
    rabit$logP <- -log10(rabit$p.value)
    rabit <- rabit[1:n, colnames(rabit) %in%  c("TF","logP")]
    rabit$TF <- factor(rabit$TF,levels = rev(rabit$TF))
    p <- ggplot(rabit, aes(x=TF,y=logP)) + geom_histogram(stat = "identity",fill="black")+theme_bw()
    p <- p + scale_y_continuous(position = "right")
    p <- p + xlab("")+ylab("-Log10(p-value)")
    p <- p+coord_flip()
    p <- p + theme(axis.text.x = element_text(size =10))
    p <- p + theme(axis.text.y = element_text(size = 10))
    p <- p + theme(axis.title.y = element_text(size = 15,face = "bold"))
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"))
    ggsave(paste(output.path,paste(str_replace(i,".txt",""),"TF.pdf",sep="_"),sep = "/"),p,height = 5,width=6)
}


#setwd("/Users/junge/Desktop/pepline_test/")
#data_path  <- "Result/SPOP/Rabit"  # join(config["output"],"Rabit")
#output_path  <- "Result/SPOP/Rabit"  # join(config["output"],"Rabit")
#data_path  <- paste(data_path,"output",sep="/")
output.path  <- paste(output_path,"view",sep="/")
dir.create(output.path,recursive = T)
for(i in list.files(data_path)){
  rabitres.path <- paste(data_path,i,sep = "/") 
  RabitView(rabitres.path=rabitres.path,output.path=output.path)
}





