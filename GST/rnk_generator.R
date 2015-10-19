#--- Authorship
## AUTHOR : Perales-Paton, Javier
## E-MAIL : jperales@cnio.es
## DATE : Oct2015
#--- About this
## TITLE : Rank (.rnk) file generator for Count-based data from RNA-seq
## DESCRIPTION : 
#   Gene set testing for RNA-seq data is not a straightforward analysis. 
# This is a R function to create a .RNK ready to be used in GSEA Java software (Broad Institute).
## DEPENDENCIES :
#   - edgeR
## It is a free piece of code. Please, let me know improvements and ideas.


# mk_rnk :: read a data.frame and make a ranked file but randomizing duplicated values.
mk_rnk <- function(df,col.ids=NULL) {
  
  # Sanity check
  if(!"data.frame" %in% is(df))
    stop("ERROR: DF arg must be a data.frame");
  
  if(is.null(col.ids) & dim(df)[2]!=2) {
    stop("ERROR: df arg must be a data.frame with two columns when col.ids is null");
    
  } 
  
  if(!is.null(col.ids) & length(col.ids)!=2)
    stop("ERROR: col.ids must be NULL or a two-length vector");
  
  if(!is.null(col.ids) & is.character(col.ids)) {
    if(!all(col.ids%in%colnames(df))) {
      stop("ERROR : colnames in col.ids were not recognized.");
    }
  }
  #
  
  if(!is.null(col.ids))
    df <- df[,col.ids]
  
  if(!any(c("factor","character")%in%is(df[,1])))
      stop("ERROR : first column must be a factor or character vector");
    
  if(!"numeric"%in%is(df[,2]))
      stop("ERROR : Second column must be a numeric vector");
  
  
  
  df.rnk <- data.frame(Symbol=toupper(df[,1]), Rank=df[,2])
  rownames(df.rnk) <- rownames(df[,1])
  
  df.rnk <- df.rnk[order(df.rnk$Rank,decreasing = F),]
  rownames(df.rnk) <- 1:nrow(df.rnk)
  
  # Randomized duplicated values
  dummy.ids <- 1:nrow(df)
  names(dummy.ids) <- 1:nrow(df)
  
  
  duplicated.values <- df.rnk$Rank[duplicated(df.rnk$Rank)]
  for(x in duplicated.values) {
    dup.ids <- which(df.rnk$Rank==x)
    dummy.ids[dup.ids] <- sample(dup.ids)
  }

  # Apply the sampling to avoid AAA,AAB,ABB order of genes 
  # which are correlated with family members (e.g. ABC proteins,OLFRXX)
  df.rnk <- df.rnk[dummy.ids,]
  
  return(df.rnk);
}

write.rnk <- function(df.rnk,output) {
  # Sanity check
  if(!"data.frame" %in% is(df.rnk))
    stop("ERROR: df.rnk arg must be a data.frame");
  
  if(dim(df.rnk)[2]!=2)
    stop("ERROR: df.rnk arg must be a data.frame with two columns");
  
  if(!any(c("factor","character")%in%is(df.rnk[,1])))
    stop("ERROR : first column must be a factor or character vector");
  
  if(!"numeric"%in%is(df.rnk[,2]))
    stop("ERROR : Second column must be a numeric vector");
  #
  
  write.table(df.rnk,file = output,sep="\t",quote = FALSE, row.names=FALSE, col.names=FALSE)
}

# PENDING TASKS:
# mk_rnk.DGEList <- function()
# mk_rnk.DGElrt <- function()