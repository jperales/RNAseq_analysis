#--- Authorship
## AUTHOR : Perales-Paton, Javier
## E-MAIL : jperales@cnio.es
## DATE : Oct2015
#--- About this
## TITLE : Cross-species contamination plot
## DESCRIPTION : 
# Create a bar plot summarizing the distribution of reads mapping on CDS, UTR, Intros, TES and TSS.
# Input: a vector of elements name=item as Sample_Name=filename from RSeQC read_distribution.py script.

## DEPENDENCIES :
#   - ggplot2
## It is a free piece of code. Please, let me know improvements and ideas.

#--- Example of input:
# files <- c(WELL_56="./Well_56_DAPI_Neg/Well_56_DAPI_Neg.read_distribution.txt",
#            WELL_64="./Well_64_DAPI_Neg/Well_64_DAPI_Neg.read_distribution.txt",
#            WELL_88="./Well_88_DAPI_Neg/Well_88_DAPI_Neg.read_distribution.txt",
#            WELL_32="./Well_32_DAPI_Neg/Well_32_DAPI_Neg.read_distribution.txt",
#            WELL_96="./Well_96_DAPI_Neg/Well_96_DAPI_Neg.read_distribution.txt",
#            WELL_27="./Well_27_DAPI_pos/Well_27_DAPI_pos.read_distribution.txt",
#            WELL_72="./Well_72_DAPI_Neg/Well_72_DAPI_Neg.read_distribution.txt",
#            WELL_24="./Well_24_DAPI_Neg/Well_24_DAPI_Neg.read_distribution.txt",
#            WELL_19="./Well_19_NoCell/Well_19_NoCell.read_distribution.txt",
#            WELL_15="./Well_15_DAPI_pos/Well_15_DAPI_pos.read_distribution.txt",
#            WELL_48="./Well_48_DAPI_Neg/Well_48_DAPI_Neg.read_distribution.txt")
# 
# DF <- read_distrib.df(files=files)
# 
# read_distrib.bar(data = DF)

#---Code

read_distrib.df <- function(files) {
  # Sanity check
  if(!"character"%in%is(files))
    stop("'files' arg must be a vector of characters");
  if(!all(file.exists(files)))
    stop("'files' does not exist");
  if(any(sapply(files,function(x) length(readLines(x)))!=16))
    stop("'files' have bad format");
  if(!all(sapply(files,function(x) grep("^=+$",readLines(x))==c(4,16))))
    stop("'files' have bad format (#2)");
  #
  
  # Fill data.frame
  DF <- vector()
  for(sname in names(files)) {
    fl.ln <- gsub("\\s+"," ",readLines(files[sname]))
    
    df <- data.frame(sname=rep(sname,10),
                   group=sapply(fl.ln[6:15],function(x) strsplit(x,split=" ")[[1]][1]),
                   count=sapply(fl.ln[6:15],function(x) as.numeric(strsplit(x,split=" ")[[1]][3])),
                   perc=100*sapply(fl.ln[6:15],function(x) as.numeric(strsplit(x,split=" ")[[1]][3]))/
                     sum(sapply(fl.ln[6:15],function(x) as.numeric(strsplit(x,split=" ")[[1]][3]))))
    DF <- rbind(DF,df);
  }
  
  DF$group <- relevel(DF$group,ref="CDS_Exons");
  
  return(DF);
}

read_distrib.bar <- function(data) {
  # Sanity check
  if(dim(data)[2]!=4)
    stop("'data' has an invalid format: bad number of columns");
  
  if(!all(colnames(data)==c("sname","group","count","perc")))
    stop("'data' has an invalid format: bad columnames");
  
  levels.format <- c("CDS_Exons","3'UTR_Exons","5'UTR_Exons","Introns",
                     "TES_down_10kb","TES_down_1kb","TES_down_5kb",
                      "TSS_up_10kb","TSS_up_1kb","TSS_up_5kb")
  if(!all(levels(data$group)==levels.format))
    stop("'data' has an invalid format: different levels")
  #
  
  require(ggplot2)
  gg <- ggplot(data,aes(sname,fill=group,y = perc)) + geom_bar(stat="identity") + theme_bw() + 
    scale_fill_manual(values = c("#00441b","#253494","#377eb8","#e41a1c","#4daf4a","#984ea3",
                                 "#ff7f00","#ffff33","#a65628","#f781bf")) +
    ylab("%Mapped Reads") + xlab("Samples") +
    theme(axis.title.y = element_text(size=14),axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=14), 
          axis.text.x = element_text(angle = 45, hjust = 1,size = 12))
  return(gg)
}

