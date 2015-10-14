#--- Authorship
## AUTHOR : Perales-Paton, Javier
## E-MAIL : jperales@cnio.es
## DATE : Oct2015
#--- About this
## TITLE : Batch effect removal
## DESCRIPTION : 
#   Correct batch effect in Expression profiles from experiments with batch effects using linear models.
## DEPENDENCIES :
#   - lattice
## It is a free piece of code. Please, let me know improvements and ideas.

#---Code

# cpm_woBatchEffect : it uses a lm function to remove batch effect and obtain a matrix of CPM per gene
# Inspirated by https://support.bioconductor.org/p/54594/
# My experience: the obtained fixed CPM and the original ones are very close (low batch effect).

# @ usage :: cpm_woBath
# @ DGEList : a DGEList from edgeR. So it is accessible the DGEList$samples$group to create the virtual design matrix.
# @ A vector grouping by batcheffect (e.g. 111222).
# @ More info: See ?removeBatchEffect
cpm_woBatchEffect <- function(DGEList,batch.vector) {
  
  if(is(DGEList)[1] != "DGEList") 
    stop("ERROR : it only works with DGEList from now because use normalized.lib.sizes factors");
  
  logCPM <- cpm(DGEList, log=TRUE, normalized.lib.sizes = T)
  
  # Design
  des <- model.matrix(~DGEList$samples$group)
  
  logCPM_batchremoved <- removeBatchEffect(logCPM,batch=batch.vector,design=des)
  CPM_batchremoved <- 2^logCPM_batchremoved
  return(CPM_batchremoved)
}