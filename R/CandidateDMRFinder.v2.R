#' CandidateDMRFinder.v2
#' 
#'    
#'    This function identifies candidate/putative differentially methylated loci (DML)
#'    based on the procedure described in Koestler et al., (2016).  Breifly, a series
#'    of two-sample t-tests are fit to the J CpGs contained in the referenceBetas 
#'    object and used to compare the mean methylation beta-values between each of the
#'    K cell type against the mean methylation beta-values computed across the remaining 
#'    K - 1 cell types.  Putative DMLs are identified by first rank ordering CpGs by their 
#'    t-statistics, then taking the top M DMLs with the smallest and largest t-statistics 
#'    for each of the K comparisons.  
#'
#' 
#'
#' 
#' @param 
#'   cellTypes:       A vector of length K that contains the cell type names.  For example,
#'					 c("CD4T", "CD8T", "NK", "Bcell", "Mono", "Gran").
#' @param                   
#'	referenceBetas:  A J x N matrix of cell-specific methylation beta-values; J represents 
#'                   the number of CpGs (i.e., ~ 450,000 for the Illumina HumanMethylation450
#'                    array) and N represents the number of samples for which cell-specific
#'                    methylation signatures are available.
#'
#'@param
#'	referenceCovars: A N x P data.frame of meta data aross the N samples.  The rows of this
#'                    object MUST be in the same order as the columns of referenceBetas.
#'                    Further, there must be a column called "CellType" (case sensitive), 
#'                    that indicates the cell-type identity for each of the N samples.
#'                    The nomenclature used to indicate cell identity across the N samples
#'                    should follow the nomenclature used for cellTypes (see above).
#'@param
#'	M:				 The number of candidate DMLs with the smallest and largest t-statistic
#'                    to return for each comparison.  Defaults to M = 150 as in Koestler et al.,
#'                    (2016)
#'@param
#'	equal.variance:  Should a t-test assuming equal variances be fit.  Defaults to FALSE, 
#'                    an unequal variance t-test.
#'@return
#'            A list containing two objects: (1) candidateSet - a vector containing the names 
#'            of the R candidate DMLs identified from the analysis and (2) coefEsts - A R x K
#'            matrix of the within-cell type mean methylation beta values across the R 
#'            identified candidate DMLs.
#' @import 	genefilter   
#' @export    


CandidateDMRFinder.v2 =function(cellTypes, referenceBetas, referenceCovars, M = 150, equal.variance = F){
    
    require(genefilter)
    
    p = referenceBetas
    pd = referenceCovars
    K = length(cellTypes)
    
    if(sum(cellTypes %in% pd$CellType)!= K) {
        stop("cell type names in target covariate data are not consistent with the nomenclature used in cellTypes")
    }
    
    splitit <- function(x) {
        split(seq(along = x), x)
    }
    
    keep <- which(pd$CellType %in% cellTypes)
    pd <- pd[keep, ]
    p <- p[, keep]
    
    tIndexes <- splitit(pd$CellType)
    tstatList1 <- lapply(tIndexes, function(i) {
        x1 <- i
        x2 <- c(1:ncol(p))[-x1]
        return(fastT(p, x1, x2, var.equal = equal.variance))
    })
    
    probeList1 <- lapply(tstatList1, function(x) {
        yUp <- rownames(p)[order(x$z, decreasing = TRUE)]
        yDown <- rownames(p)[order(x$z, decreasing = FALSE)]
        c(yUp[1:M], yDown[1:M])
    })
    
    candidateSet <- unique(unlist(probeList1))
    p <- p[candidateSet, ]
    
    coefEsts = matrix(NA, nrow = length(candidateSet), ncol = K)
    rownames(coefEsts) = candidateSet
    colnames(coefEsts) = cellTypes
    for(k in 1:K) {
        ind = which(pd$CellType %in% cellTypes[k])
        coefEsts[,k] = apply(p[,ind], 1, mean, na.rm = T)
    }
    
    tmp = list(candidateSet, coefEsts)
    names(tmp) = c("candidateSet", "coefEsts")
    tmp
}