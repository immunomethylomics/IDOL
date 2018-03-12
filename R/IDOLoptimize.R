#' IDOLoptimize
#' 
#' 
#' 
#'@description   This function identifies optimal DML/DMR libraries for cell mixture deconvolution using 
#'    the procedure described in Koestler et al., (2016).  
#'
#'
#'@param   candDMRFinderObject:    List object returned from the CandidateDMRFinder.v2 function 
#'                    
#'@param	trainingBetas:  		A J x N matrix of beta-values where methylation was profiled 
#'                           in a heterogenous tissue type (i.e., WB).  Here, J indicates
#'                           the number of CpGs and N, the number of samples
#'
#'@param	trainingCovariates: 	A N x P data.frame of meta data aross the N samples.  Contained
#'                           within the meta data must be the observed cell fractions (i.e. FACS or
#'                           otherwise) for the K cell types indicated in the list object returned 
#'                           from the CandidateDMRFinder.v2 function.  The naming of cell types
#'                           should be consistent between these two objects as well.
#'                           
#'@param	libSize:				Size of the optimized IDOL library.  Defaults to 300.
#'
#'@param	maxIt:  				Maximum number of iterations for the IDOL algorithm. Defaults 
#'                           to 500.
#'
#'@param   numCores:               Number of processing cores to use (see R package DoParallel).  
#'                           Defaults to 4.
#'
#'@return
#'    A list containing six objects: (1) "IDOL Optimized Library" - a vector containing the names 
#'            of the CpGs in the identified IDOL optimized library (2) "IDOL Optimized CoefEsts" - 
#'            matrix of the within-cell type mean methylation beta values across the CpGs in the
#'            optimal library (3) "RMSE" average root-mean squared error calculated across each 
#'            iteration of IDOL (4) "R2" average R2 (coefficent of determinatino) calculated across each 
#'            iteration of IDOL (5) "Number of iterations" how many iterations of IDOL were used
#'            (6) "Library Size" libsize above.
#' @import quadprog  
#' @import doParallel     
#'@export 


IDOLoptimize = function(candDMRFinderObject, trainingBetas, trainingCovariates, 
                        libSize = 300, maxIt = 500, numCores = 4) {
    
    # Load necessary packages
    require(quadprog)
    require(doParallel)
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    
    # Define relevant functions
    expit = function(w) exp(w)/(1 + exp(w))
    logit = function(w) log(w) - log(1-w)
    
    R2compute = function(obs, pred) {
        r2 = NULL
        for(i in 1:dim(obs)[2]) {
            y = obs[,i]
            x = pred[,i]
            r2[i] = summary(lm(y~x))[[8]]
        }
        sum(r2)/dim(obs)[2]
    }
    
    # Polar coordinates function 
    polar = function(x, y, scale = 1) {
        r = sqrt(x^2 + y^2)
        theta = atan2(y, x)
        r*cos(theta - (scale*pi/4))
    }
    
    # Define relevant parameters
    trainingProbes1 = candDMRFinderObject$candidateSet
    coefEsts = candDMRFinderObject$coefEsts
    P = length(trainingProbes1)
    ProbVector = rep(1/P, P)
    V = libSize  # number of CpGs to select
    B = maxIt # number of interations
    
    R2.null = 0
    R2Vals = NULL
    RMSE.null = 10000
    RMSEVals = NULL
    
    # do some cross-checks before beginning algorithm
    cellTypes = colnames(coefEsts)
    K = length(cellTypes)
    if(sum(cellTypes %in% colnames(trainingCovariates)) != K) {
        stop("cell type names in target covariate data are not consistent with cell types to be deconvoluted")
    }
    
    # Perform the idol maximization algorithm
    for(i in 1:B) {
        Probes = sample(1:P, V, prob = ProbVector)
        CpGNames = trainingProbes1[Probes]
        
        Beta = coefEsts[CpGNames, ]
        Lwbc = diag(ncol(Beta))
        
        ctpred = data.frame(minfi:::projectCellType(trainingBetas[CpGNames,], Beta))
        omega.tilde = 100*ctpred
        omega.obs = trainingCovariates[,cellTypes]
        
        # compute RMSE based on all probes
        # Whole blood
        diff = as.matrix(omega.tilde - omega.obs)
        RMSE = sqrt(sum(diff^2)/K)
        
        # compute R2 based on all probes
        R2 = R2compute(omega.obs, omega.tilde)
        
        # compute the leave-one-out R-squared values for each of the probes
        Perform.q = foreach(j = 1:length(CpGNames)) %dopar% {
            
            Beta.q = Beta[CpGNames[-j],]
            
            ctpred.q = data.frame(minfi:::projectCellType(trainingBetas[CpGNames[-j],], Beta.q, Lwbc))
            omega.tilde.q = 100*ctpred.q
            
            # compute RMSE based on all probes
            # Whole blood
            diff.q = as.matrix(omega.tilde.q - omega.obs)
            RMSE.q = sqrt(sum(diff.q^2)/K)
            
            # compute R2 based on all probes
            R2.q = R2compute(omega.obs, omega.tilde.q)
            cbind(R2.q, RMSE.q)
        }
        
        R2.q = unlist(Perform.q)[seq(from = 1, to = (2*V-1), by = 2)]
        RMSE.q = unlist(Perform.q)[seq(from = 2, to = 2*V, by = 2)]
        
        # compute relevant values to modify the probabilities that a probe is selected in subsequent iterations
        rmse.dq = (RMSE - RMSE.q)*(-1)
        norm.rmse = (rmse.dq)/sd(rmse.dq)
        
        r2.dq = (R2 - R2.q)
        norm.r2 = (r2.dq)/sd(r2.dq)
        
        p1 = polar(norm.rmse, norm.r2)
        
        for (j in 1:length(Probes)) {
            p0 = ProbVector[[Probes[j]]]
            ProbVector[[Probes[j]]] = expit(p1[j])*p0 + p0/2
        }
        
        # rescale the ProbMat
        ProbVector = ProbVector/sum(ProbVector)
        
        # save the optimal library for later use
        if(RMSE <= RMSE.null & R2 >= R2.null) {
            RMSE.null = RMSE
            R2.null = R2
            print(paste("Iteration: ", i, " RMSE=", round(RMSE, 3), "; R2=", round(R2, 3), sep = ""))
            
            IDOL.optim.DMRs = CpGNames
            IDOL.optim.coefEsts = coefEsts[CpGNames,]
            
            save(IDOL.optim.DMRs, IDOL.optim.coefEsts, file = paste("IDOL optimized DMR library_", V, ".RData", sep = ""))
        }
        RMSEVals[i] = RMSE
        R2Vals[i] = R2
    }
    stopCluster(cl)
    IDOLObjects = list(IDOL.optim.DMRs, IDOL.optim.coefEsts, RMSEVals, R2Vals, B, V)
    names(IDOLObjects) = c("IDOL Optimized Library", "IDOL Optimized CoefEsts", 
                           "RMSE", "R2", "Number of Iterations", "LibrarySize")
    
    print(paste("The Average RMSE = ", round(RMSE.null, 3), " and R2 = ", round(R2.null, 3), 
                " for the IDOL Optimized Library", sep = ""))
    
    return(IDOLObjects)
}





