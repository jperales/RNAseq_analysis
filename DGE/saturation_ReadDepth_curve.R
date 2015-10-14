
library("Rsamtools")

which <- RangesList(seq1=IRanges(1000, 2000), seq2=IRanges(c(100,1000),c(1000,2000)))
what <- c("rname","strand","pos","qwidth","seq")

param <- ScanBamParam(which=which, what=what)
bamFile <-system.file("extdata", "ex1.bam", package="Rsamtools")
bam <- scanBam(bamFile, param=param)


.unlist <- function (x) {
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)) {
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
}

bam <- unname(bam)
elts <- setNames(bamWhat(param),bamWhat(param))
lst <- lapply(elts, function(elt) .unlist(lapply(bam, "[[", elt)))
head(do.call("DataFrame", lst))
