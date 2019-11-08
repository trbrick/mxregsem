#
#   Copyright 2007-2019 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#' External search for an optimal penalty value
#'
#' Perform a grid search over a space of regularization penalties.
#' This happens in the front-end, returning to R at each cycle.
#' Currently only guaranteed to work for single-parameter models.
#' Expect many changes in the interface while this is developed further.
#' **WARNING: Dangerously assumes that fitfunction is the -2 log likelihood**
#' 
#' The EBIC computation uses the computation from the previous regsem.
#' It is not clear why this is the case.
#'
#' @param model an MxRegularizedModel object to regularize
#' @param epsilon how close to zero is zero?
#' @param ... the search space for this sequence in format *param=c(...)* where c(...) is the list of possible values for parameter param.
#' @param approach what fit function to use?  Currently only EBIC is available
#' @param ebicGamma what Gamma value to use for EBIC?  Must be between 0 and 1
#' @param returnConstrained fix all regularized-to-zero values to zero, and rerun fit without regularization
#' 
#' @details TODO: Figure out how to return all the individual runs.
#' 
#' @return A regularized MxModel; at this point an MxRegularizedModel
#' @import OpenMx
# @importFrom methods .hasSlot is slotNames
#' @importFrom stats sd
#' @export
#' 
mxPenaltySearchExternal <- function(model, search_space=NULL, epsilon=1e-6, ..., approach="EBIC", ebicGamma=.5, returnConstrained=FALSE, verbose=TRUE) {
  
  # Check for regularization
  if(!is(model, "MxRegularizedModel") || length(model@regularizations) < 1) {
    warning("Not a regularized model.")
    return(mxRun(model))
  }
  
  # Find all the hyperparameters in each penalty object
  penalties <- model@regularizations
  for(penalty in penalties) {
    penaltySearch <- NULL
    thisName <- penalty@name
    penMat <- paste0(thisName, "_Penalty")
    for(aVal in names(penalty@hyperparameters)) {
      if(is.null(search_space) || aVal %in% names(search_space)) {
        if(is.null(search_space)) {
          thisParam <- data.frame(penalty@hpranges[[aVal]])
        } else {
          thisParam <- data.frame(search_space[[aVal]])
        }
        names(thisParam) <- c(paste0(thisName, "_", aVal))
        if(is.null(penaltySearch)) {
          penaltySearch <- thisParam
        } else {
          penaltySearch <- suppressWarnings(cbind(penaltySearch[rep(1:nrow(penaltySearch), length(thisParam)), 1:ncol(penaltySearch), drop=FALSE], thisParam[rep(1:nrow(thisParam), each=nrow(penaltySearch)), 1:ncol(thisParam),drop=FALSE])) # Suppress to avoid row-name-dropping warnings
        }
      }
    }
  }

  getEBIC <- function(idx, model, search_space, ..., gamma=ebicGamma, pb=NULL) {
    amodel <- omxSetParameters(model, labels=names(search_space), values=as.numeric(search_space[idx,]))
    amodel <- mxOption(amodel, "Calculate Hessian", "No")
    amodel <- mxOption(amodel, "Standard Errors"  , "No")
    tFit <- mxRun(amodel, silent=TRUE, unsafe=TRUE) #, beginMessage=FALSE, suppressWarnings=TRUE)
    tSummary <- summary(tFit, epsilon=epsilon, verbose=FALSE, refModels=NULL, computeRefs=FALSE)
    EP <- tSummary$estimatedParameters
    DF <- tSummary$degreesOfFreedom
    sub1 <- tFit@submodels[[1]]
    sub1name <- tFit@submodels[[1]]$name
    ll <- mxEvalByName(paste0(sub1name, ".fitfunction"), tFit)[1,1]  # Danger: assumes fit function is -2LL
    N <- tSummary$numObs
    p <- length(sub1$manifestVars)
    nfac <- length(sub1$latentVars)
    
    # EBIC -- Need a better function.
    if(nfac < 1) {
      nfac <- 1
    }
    if(model$data$type == "raw") {
      EBIC <- (ll)/((N)/(p+1)) + (log(N) * EP + 2*EP * gamma * log(p+nfac))
    } else {
      EBIC <- (ll) + log(N) * EP + 2*EP * gamma * log(p + nfac)
    }
    params <- omxGetParameters(tFit)
    # Progress!
    if(!is.null(pb)) {
      setTxtProgressBar(pb, idx)
    }
    
    # output. TODO: Something broken here?
    return(c(EBIC=EBIC, ll=ll, EP=EP, DF=DF, gamma=gamma, nfac=nfac, p=p, 
             omxGetParameters(tFit, labels=names(search_space) )))
  }
  
  # Start at good starting values.
  model <- mxTryHard(model, extraTries = 100)
  
  # Progress bar setup.
  pb=NULL
  if(verbose==TRUE){
    pb <- txtProgressBar(min = 0, max = nrow(penaltySearch), style = 3)
  }
  
  outList <- sapply(1:nrow(penaltySearch), getEBIC, search_space=penaltySearch, epsilon=epsilon, model=model, gamma=ebicGamma, pb=pb)
  
  # Name preservation is apparently difficult.
  tNames <- rownames(outList)
  outList <- data.frame(t(outList))
  names(outList) <- tNames
  
  if(returnConstrained) {
    minmod <- unlist(outList[which.min(outList$EBIC),])
    fixed <- abs(minmod) < epsilon
    minmod[fixed] <- 0
    names(fixed) <- names(minmod)
    tries <- outList
    tParms <- omxGetParameters(model, free=TRUE)
    combinedNames <- intersect(names(minmod), names(tParms))
    tmodel <- omxSetParameters(model, labels=combinedNames, 
                               values=minmod[combinedNames], 
                               free=!fixed[combinedNames])
    tmodel <- OpenMx::mxModel(tmodel, mxComputeDefault())
    outList <- mxTryHard(tmodel)
    otherNames <- setdiff(names(minmod), names(tParms))
    message("Final solution has params ", paste(names(omxGetParameters(outList)), "=", omxGetParameters(outList), collapse=", "), " and a likelihood of ", paste0(summary(outList)$Minus2LogLikelihood), " at reg levels ", paste0(otherNames, ":", minmod[otherNames], collapse=", "), "and an N of ", outList$data$numObs)
    outList@output <-c(outList@output, list("regularizationSearchOutcomes"=tries))
  } else {
    minmod <- unlist(outList[which.min(outList$EBIC),])
    tries <- outList
    tRegs <- names(penaltySearch)
    tParms <- omxGetParameters(model, free=TRUE)
    combinedNames <- intersect(names(minmod), names(tParms))
    combinedNames <- setdiff(combinedNames, tRegs)
    # Filter unnamed, which somehow choke everything.
    combinedNames <- setdiff(combinedNames, grep("[[]", combinedNames, value = TRUE))
    tmodel <- omxSetParameters(model, labels=combinedNames, 
                               values=minmod[combinedNames], 
                               free=TRUE)
    tmodel <- omxSetParameters(tmodel, labels=tRegs, 
                               values=minmod[tRegs], 
                               free=FALSE)
    outList <- mxTryHard(tmodel)
    otherNames <- setdiff(names(minmod), names(tParms))
    message("Final solution has params ", paste(names(omxGetParameters(outList)), "=", omxGetParameters(outList), collapse=", "), " and a likelihood of ", paste0(summary(outList)$Minus2LogLikelihood), " at reg levels ", paste0(otherNames, ":", minmod[otherNames], collapse=", "), "and an N of ", outList$data$numObs)
    outList@output <-c(outList@output, list("regularizationSearchOutcomes"=tries))
  }
  
  # Trying to reintegrate into an openMx model output.
  # minmod <- unlist(outList[which.min(outList$EBIC),])
  # tParms <- omxGetParameters(model, free=TRUE)
  # combinedNames <- intersect(names(minmod), names(tParms))
  # model <- omxSetParameters(model, labels=combinedNames, values=minmod[combinedNames])
  # model <- mxModel(model, mxComputeSequence(list(
  #                           mxComputeOnce(paste(sub1name, 
  #                                               c('expectation',
  #                                                 'fitfunction')
  #                           mxComputeStandardError(),)))

  return(outList)
} 




#' Search for an optimal penalty value  **Experimental**
#'
#' Perform a grid search over a space of regularization penalties.
#' Currently only guaranteed to work for single-parameter models.
#' Expect many changes in the interface while this is developed further.
#'
#' @param model an MxRegularizedModel object to regularize
#' @param ... the search space for this sequence in format *param=c(...)* where c(...) is the list of possible values for parameter param.
#' @return A regularized MxModel; at this point an MxRegularizedModel
#' @import OpenMx
# @importFrom methods .hasSlot is slotNames
#' @importFrom stats sd
#' @export
#' 
mxPenaltySearch <- function(model, search_space=NULL, epsilon=1e-6, ...) {
  stop("Error NYI: Internal search is not yet available. Use mxPenaltySearchExternal() for now.")
  # Check for regularization
  if(!is(model, "MxRegularizedModel")) {
    stop("Not a regularized model.")
  }
  
  # EBIC computation
  UncorrectedDF <- mxMatrix()
  EBIC_string <- paste0(model$submodels[[1]]$name + ".fitfunction + log(", model$submodels[[1]]$data$numObs, ") %*% ")
  model@algebras[["EBIC_EVAL"]] <- mxAlgebraFromString(model$submodels[[1]]$name)
  
  # Prep the base compute list
  
  
  # Find all the hyperparameters in each penalty object
  penalties <- model@regularizations
  for(penalty in penalties) {
    penaltySearch <- data.frame()
    thisName <- penalty@name
    penMat <- paste0(thisName, "_Params")
    for(aVal in names(penalty@hyperparameters)) {
      if(aVal %in% names(search_space)) {
        thisParam <- search_space[[aVal]]
        penaltySearch <- cbind(penaltySearch[rep(1:nrow(penaltySearch), length(thisParam)),], TEMPNAME=rep(thisParam, each=nrow(penaltySearch)))
        names(penaltySearch)[ncol(penaltySearch)-1] <- aVal
      }
    }
    # Add this to the compute list
  }
  
  # Complete compute list
  model <- mxModel(model, 
   mxComputeSequence(steps=list(
    mxComputeOnce('fitfunction', 'fit'),
    mxComputeNumericDeriv(parallel=FALSE, iterations=2L),
    mxComputeStandardError(),
    mxComputeHessianQuality(),
    mxComputeReportDeriv())
   )
  )
  
}
