#--- Authorship
## AUTHOR : Perales-Paton, Javier
## E-MAIL : jperales@cnio.es
## DATE : Oct2015
#--- About this
## TITLE : Cross-species contamination plot
## DESCRIPTION : 
#     Xenografts and other technical problems could arise with cross-species contamination in your raw data,
#   this script try to detect these issues.

## DEPENDENCIES :
#   - lattice
## It is a free piece of code. Please, let me know improvements and ideas.

#---Code

# crossSps.readout :: Get a matrix of values to plot by parsing fastq_screen output files
#@ usage : files : char vector where name element is the sample name, and the element is the file
crossSps.readout <- function(files) {
  if(!all(c("vector","character") %in% is(files)))
    stop(paste0("'files' argument must be a character vector indicating path to files."))
  
  # Columns of interest in fastq_screen output files
  cols.int <- c(5,7,9,11) # After removing the first element, these colIds are related to : 
  # 5: %One_hit_one_library 
  # 7: %Multiple_hits_one_library
  # 9: %One_hit_multiple_libraries
  # 11: %Multiple_hits_multiple_libraries
  
  
  mt <- vector()
  for(i in 1:length(files)) {
    if(file.exists(files[i])) {
        lns <- readLines(files[[i]]) # It is a messy table, so read lines and parse it
        # Sanity check : check out fastq_screen format
        if(!grepl("#Fastq_screen",lns[1]))
          stop("Bad format in first line.");
        
        if(lns[grep("%Hit_no_libraries",lns)-1]!="") 
          stop("Bad format : end not recognized");
        #
        
        # Parsing
        Nlibs <- grep("%Hit_no_libraries",lns)-2; # How many libs you tried out given your configuration
        
        mt.tmp <- sapply(lns[3:Nlibs],function(x) as.numeric(strsplit(x,split="\t")[[1]][-1]))  # It works out well
        colnames(mt.tmp) <- sapply(colnames(mt.tmp),function(x) strsplit(x,split = "\t")[[1]][1])
        mt.tmp <- mt.tmp[cols.int,]
        mt.tmp <- cbind(mt.tmp,
                        No_hit=c(as.numeric(strsplit(grep("%Hit_no_libraries",lns,value=TRUE),split=": ")[[1]][2]),
                                 rep(0,nrow(mt.tmp)-1)))
        rownames(mt.tmp) <- gsub("%","",strsplit(lns[[2]],split="\t")[[1]][-1][cols.int])
        mt2.tmp <- t(mt.tmp)
        mt3.tmp <- as.data.frame(mt2.tmp)
        mt.sample <- cbind(SName=rep(names(files)[i],nrow(mt3.tmp)),Org=rownames(mt3.tmp),mt3.tmp)
        mt <- rbind(mt,mt.sample)
    } else {
      stop(paste0("Critical error : ",files[[i]]," File does _NOT_ exist. Aborting..."))
    }
  }
  # Reorder by original order
  mt$Org <- factor(mt$Org,levels(mt$Org)[match(colnames(mt.tmp),levels(mt$Org))])
  
  return(mt)
}


crossSps.plot <- function(crossSps.matrix,nrow.grid=4,ncol.grid=4) {
  # Function to check if it is integer
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # Args sanity check
  mandatory.cols <- c("SName","Org",
                      "One_hit_one_library","Multiple_hits_one_library",
                      "One_hit_multiple_libraries","Multiple_hits_multiple_libraries")
  if(!"data.frame"%in%is(crossSps.matrix) | ncol(crossSps.matrix)!=6 |
     !all(mandatory.cols%in%colnames(crossSps.matrix)))
     stop("Bad data.frame in 'crossSps.matrix");
  
  if(!(is.numeric(nrow.grid) & length(nrow.grid)==1))
    stop("'nrow.grid' must be numeric");
  if(!(is.numeric(ncol.grid) & length(ncol.grid)==1))
    stop("'ncol.grid' must be numeric");
  #
  if(!is.wholenumber(nrow.grid) | !is.wholenumber(ncol.grid))
    stop("'ncol.grid' and 'nrow.grid' must be")
  
  # Libs
  require(lattice)
  # require(RColorBrewer)
  
  # myColours.blue <- brewer.pal(6,"Blues")
  myColours.blue <- C("#EFF3FF","#C6DBEF","#9ECAE1","#6BAED6","#3182BD","#08519C")
  # myColours.red <- brewer.pal(6,"Reds")
  myColours.red <- c("#FEE5D9","#FCBBA1","#FC9272","#FB6A4A","#DE2D26","#A50F15")
  
  # Lattice settings:
  my.settings <- list(
    superpose.polygon=list(col=c(myColours.blue[c(5,2)],myColours.red[c(5,2)]), border="black"),
    strip.background=list(col=myColours.blue[6]),
    # layout.widths = list(right.padding = 5)
    strip.border=list(col="black")
  )
  
  # Plot
  barchart(One_hit_one_library+Multiple_hits_one_library+One_hit_multiple_libraries+Multiple_hits_multiple_libraries~Org|SName,
           data=crossSps.matrix,stack=TRUE,layout=c(ncol.grid, nrow.grid),
           scales=list(y=list(relation="same",rot=0,cex=1,alternating=3,tck=0.5,limits=c(0,100)),
                       x=list(rot=45,cex=1)),
           ylab=list(label="%Mapped Reads",cex=1.5),
           par.settings= my.settings,
           par.strip.text=list(col="white", font=1,cex=1.5),
           between = list(x = 0.5),
           # par.settings = list(layout.widths = list(right.padding = 5)),
           auto.key=list(space="bottom", columns=4, points=FALSE, rectangles=TRUE,cex.title=0.8))
  
}

# crossSps :: Execute the two functions to get the whole picture
#@ usage :: crossSps.fig(fls,nrow.grid,ncol.grid)

crossSps <- function(files,nrow.grid=4,ncol.grid=4) {
  mt <- crossSps.readout(files = files);
  crossSps.plot(crossSps.matrix = mt,nrow.grid=nrow.grid,ncol.grid=ncol.grid);
}