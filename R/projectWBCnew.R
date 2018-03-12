#############################################################################################
#' projectWBCnew
#' 
#' 
#' 
#'@description    This function predicts the underlying cellular composition of heterogeneous tissue 
#'    types (i.e., WB) using the constrained projection procedure described by Houseman et al., 
#'    (2012).  
#'
#'
#'      
#'@param	Y:  			 A J x N matrix of methylation beta-values collected from mixed/
#'                    heterogeneous biospecimen (i.e., WB).  Target set.
#'
#' @param 	coefWBC:         A J x K projection matrix;, i.e., within-cell type mean methylation 
#'                    matrix across J DMLs and K many cell types
#'
#'@param	contrastWBC:	 Contrast for cell composition predictions.  The user needn't modify
#'                    this 
#'
#'	nonnegative:     Should cell predictions be nonnegative.  Defaults to TRUE
#'
#'@return 
#'    A N x K matrix of cell proportion estimates across the K cell types for each 
#'            of the N subjects contained in the Target Set.
#' @import 	quadprog    
#'@export    

projectWBCnew = function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE){ 
    
    if(is.null(contrastWBC)) Xmat = coefWBC
    else Xmat = coefWBC %*% t(contrastWBC) 
    
    nCol = dim(Xmat)[2]
    nSubj = dim(Y)[2]
    
    mixCoef = matrix(0, nSubj, nCol)
    rownames(mixCoef) = colnames(Y)
    colnames(mixCoef) = colnames(Xmat)
    
    if(nonnegative){
        library(quadprog)
        
        Amat = cbind(rep(-1,nCol), diag(nCol))
        b0vec = c(-1,rep(0,nCol))
        
        for(i in 1:nSubj){
            obs = which(!is.na(Y[,i])) 
            Dmat = t(Xmat[obs,])%*%Xmat[obs,]
            mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec, meq = 0)$sol
        }
    }
    else{
        for(i in 1:nSubj){
            obs = which(!is.na(Y[,i])) 
            Dmat = t(Xmat[obs,])%*%Xmat[obs,]
            mixCoef[i,] = solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
        }
    }
    
    return(mixCoef)
}