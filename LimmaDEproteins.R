#' Differential protein abundance analysis 
#' @name LimmaDEproteins
#' @param exprSet Expression data after log2 transform
#' @param tre  Mutated or Treated condition.It should be recognized in the colnames of exprSet
#' @return Limma_DE_proteins  The result of limma.
#' @author Jun Ge
#' @import limma
LimmaDEproteins <- function(exprSet=exprSet,tre="M"){
  tre_ctrl <- factor(ifelse(grepl(tre, colnames(exprSet)),"Mutated","Ctrl"),levels = c("Mutated","Ctrl"))
  # Make design matrix
  design <- model.matrix(~0+tre_ctrl) # define matrix
  rownames(design) <- colnames(exprSet)
  colnames(design) <- levels(factor(tre_ctrl))
  # Make contrast matrix
  contrast.matrix<-makeContrasts(contrast="Mutated-Ctrl",levels = design) 
  # Limma model
  fit1 <- lmFit(exprSet, design=design) 
  fit2 <- contrasts.fit(fit1,contrasts = contrast.matrix)
  fit2 <- eBayes(fit2)
  Limma_DE_proteins <- topTable(fit2, coef=1, n=Inf)
  return(Limma_DE_proteins)
}
