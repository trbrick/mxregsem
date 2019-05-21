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

#' Regularize an MxModel
#'
#' Create a regularized container object that encloses an unregularized mxModel.
#' This was the first trial; you should now just create MxModels using MxModel.
#'
#' @param model the original mxmodel to regularize
#' @param reg_params the list of parameters or algebra cells regularized by this object.  The default (NULL) regularizes on all free parameters in the model.
#' @param penalty_function the type of function to eb used for regularization.  Currently supports "lasso" or "ridge"
#' @param ... arguments to the penalty function, usually `lambda=` (for lasso or ridge) and `alpha=` (for elastic net)
#' @return A regularized MxModel; at this point a wrapper model with regularization enclosed
#' @import OpenMx
# @importFrom methods .hasSlot is slotNames
#' @importFrom stats sd
#' @export
#' 
regularizeMxModel <- function(model, reg_params=NULL, penalty_value = 50, penalty_function="lasso") {
  
   .Deprecated("mxModel", package="mxregsem", msg=paste("regularizeMxModel()",
                           "is now deprecated. Use mxregsem's mxModel function",
                           "and add an mxPenalty to it. See vignette(reg_cfa2).")
                            )
  handledRegularizers <- c("lasso", "ridge")

  # Inputs:
  pfun <- pmatch(penalty_function, handledRegularizers)
  if(is.na(pfun)) {
    stop("Error NYI: ",
         penalty_function, " is not a supported penalty type yet.  ",
         "The list of currently supported reg functions is ", handledRegularizers)
  }

  # Doesn't handle def vars well yet.
  paramVals <- omxGetParameters(model, fetch="values")[reg_params]
  paramsFree <- TRUE
  paramsFree[grep("[\\[]", reg_params)] <- FALSE

  penalty <- NULL  # Scoping
  pAlg <- NULL     # Scoping
  if(pfun == 1) {
    penalty <- mxAlgebra(sum(abs(params)), name="lasso") # Included as separate algebra for clarity
    pAlg <- mxAlgebraFromString(paste0(model$name, ".fitfunction + ", penalty_value, " * lasso"),name="Fit")
  } else if(pfun == 2) {
    penalty <- mxAlgebra(t(params) %*% params, name="ridge") # Included as separate algebra for clarity
    pAlg <- mxAlgebraFromString(paste0(model$name, ".fitfunction + ", penalty_value, " * ridge"), name="Fit")
  } else {
    stop(paste("NYI: Regularization approach", penalty_function, "is not known.  No way to do arbitrary computations yet.  Sorry!"))
  }

  regModel <- mxModel(paste("regularized", model$name, sep="_"),
                      mxMatrix("Full", nrow=length(reg_params), ncol=1, values=paramVals, free=paramsFree,
                               labels=reg_params, name="params"),
                      penalty,
                      pAlg,
                      mxFitFunctionAlgebra("Fit"),
                      model, model$data
  )
  return(regModel)
}


#' getParamsInMatrix
#'
#' Get all free parameters from a given matrix
#'
#' @param model the MxModel
#' @param matrix The matrix from which to draw all associated
#' @return a vector containing the names of parameters in a given matrix in the model
#' @import OpenMx
#' @export
getParamsInMatrix <- function(model, matrix) {
  # browser()
  if(is.null(model[[matrix]]) && is(model, "MxRegularizedModel")) {
    # Regularized model passthrough
    model <- model$submodels[[1]]
  }
  newModel <- OpenMx::mxModel(model$name, model[[matrix]])
  prams <- omxGetParameters(newModel)
  return(names(prams))
}



#' Regularize an MxModel
#'
#' Create a regularized container object that encloses an unregularized mxModel.
#'
#' @param model the regularized model to summarize
#' @param ... values passed to mxSummarize
#' @param verbose gives more summary feedback
#' @param penalty_function the type of function used for regularization. Default guesses from regularization algebra name
#' @param epsilon how close to zero is close enough? Regularization limits for DF correction
#' @return Summary object summarizing the mxModel
#' @import OpenMx
summarizeRegularized <- function(regFit, ..., verbose=FALSE, penalty_function="guess", epsilon = 1e-6) {
  # default penalty_function ("guess") looks for a "lasso" object in the container model.
  dotArguments <- list(...)
  if(penalty_function == "guess") {
    penalty_function=c("lasso", "ridge")
  }
  if(length(regFit$submodels) == 1 && any(penalty_function %in% names(regFit$algebra))) {
    model <- regFit$submodels[[1]]
  } else {
    model <- regFit
  }
  if (!is.null(dotArguments[["refModels"]])) {
    refModels <- dotArguments[["refModels"]]
  } else {
    refModels <- OpenMx:::mxRefModels(regFit$submodels[[1]], run=TRUE)
  }
  satModel <- refModels[['Saturated']]
  indModel <- refModels[['Independence']]
  saturatedLikelihood <- OpenMx:::refToLikelihood(satModel)
  saturatedDoF <- OpenMx:::refToDof(satModel)
  independenceLikelihood <- OpenMx:::refToLikelihood(indModel)
  independenceDoF <- OpenMx:::refToDof(indModel)

  numObs <- dotArguments$numObs
  numStats <- dotArguments$numStats

  # Handle lasso regularization
  regOffset <- 0
  penFun <- pmatch(penalty_function, c("guess", "lasso"), nomatch = "100")
  if(!is.na(penFun) && (penFun == 2 || (penFun==1 && !is.null(regFit[['lasso']])))) {
    regs <- regFit[["params"]]
    if(is.null(regs)) {
      stop(paste("Error NYI: Not able to find regularization parameters in model", regFit$name))
    }
    regOffset <- sum(abs(regs$values) < epsilon)
  }

  useSubmodels <- dotArguments$indep
  if (is.null(useSubmodels)) { useSubmodels <- TRUE }
  retval <- list(wasRun=model@.wasRun, stale=model@.modifiedSinceRun)
  retval$parameters <- OpenMx:::parameterList(model, useSubmodels)
  if (!is.null(model@compute$steps[['ND']]) && model@compute$steps[['ND']]$checkGradient &&
      !is.null(model@compute$steps[['ND']]$output$gradient)) {
    gdetail <- model@compute$steps[['ND']]$output$gradient
    retval$seSuspect <- !gdetail[,'symmetric']
    names(retval$seSuspect) <- rownames(gdetail)
  }
  if (is(model@compute, "MxComputeBootstrap")) {
    bq <- c(.25,.75)
    if (!is.null(dotArguments[["boot.quantile"]])) {
      bq <- sort(as.numeric(dotArguments[["boot.quantile"]]))
    }
    summaryType <- 'bcbci'
    if (!is.null(dotArguments[["boot.SummaryType"]])) {
      summaryType <- dotArguments[["boot.SummaryType"]]
    }
    cb <- regFit@compute
    if (!is.null(cb@output$raw) && is.na(cb@only) && cb@output$numParam == nrow(retval$parameters)) {
      raw <- cb@output$raw
      mask <- raw[,'statusCode'] %in% cb@OK
      bootData <- raw[mask, 3:(nrow(retval$parameters)+2), drop=FALSE]
      if (sum(mask) < .95*nrow(raw)) {
        pct <- round(100*sum(mask) / nrow(raw))
        warning(paste0("Only ",pct,"% of the bootstrap replications ",
                       "converged. Accuracy is much less than the ", nrow(raw),
                       " replications requested"), call.=FALSE)
      }
      if (sum(mask) >= 3) {
        retval$bootstrapSE <- apply(bootData, 2, sd)
        retval$bootstrapQuantile <-
          OpenMx:::summarizeBootstrap(retval$parameters[, 'Estimate'], bootData, bq, summaryType)
      }
    }
  } else if (any(grep('^boot\\.', names(dotArguments)))) {
    warning("No bootstrap data found. See ?mxBootstrap")
  }
  retval$GREMLfixeff <- OpenMx:::GREMLFixEffList(model)
  retval$infoDefinite <- regFit@output$infoDefinite
  retval$conditionNumber <- regFit@output$conditionNumber
  if (length(regFit@output$gradient)) {
    agrad <- abs(regFit@output$gradient)
    retval$maxAbsGradient <- agrad[ order(-agrad)[1] ]
  }
  retval <- OpenMx:::boundsMet(model, retval)
  retval <- OpenMx:::setLikelihoods(model, saturatedLikelihood, independenceLikelihood, retval)
  retval <- OpenMx:::setNumberObservations(numObs, regFit@runstate$datalist, regFit@runstate$fitfunctions, retval)
  retval[["regularizedParameters"]] <- regOffset
  retval <- computeOptimizationStatistics(regFit, numStats, useSubmodels, saturatedDoF, independenceDoF, retval)
  retval$dataSummary <- OpenMx:::generateDataSummary(model, useSubmodels)
  retval$CI <- as.data.frame(regFit@output$confidenceIntervals)
  if (length(retval$CI) && nrow(retval$CI)) {
    retval$CI <- cbind(retval$CI, note=apply(retval$CI, 1, function(ci) {
      # This should probably take into account whether both bounds
      # were requested and consider the optimizer codes also. TODO
      if (any(is.na(ci)) || ci[1] == ci[3] || ci[1] >= ci[2] || ci[2] >= ci[3]) {
        "!!!"
      } else {
        ""
      }
    }))
  }
  retval$CIcodes <- regFit@output$confidenceIntervalCodes
  statusCode <- regFit@output$status$code
  if (!is.null(statusCode)) {
    message <- OpenMx:::optimizerMessages[[as.character(statusCode)]]
    retval[['npsolMessage']] <- message
    retval[['statusCode']] <- OpenMx:::as.statusCode(statusCode)
    retval[['maxRelativeOrdinalError']] <- regFit@output[['maxRelativeOrdinalError']]
  }
  if( .hasSlot(regFit,"compute") && length(regFit$compute$steps$CI) ){
    retval$CIdetail <- regFit$compute$steps$CI$output$detail
  }
  retval$timestamp <- regFit@output$timestamp
  retval$frontendTime <- regFit@output$frontendTime
  retval$backendTime <- regFit@output$backendTime
  retval$independentTime <- regFit@output$independentTime
  retval$wallTime <- regFit@output$wallTime
  retval$cpuTime <- regFit@output$cpuTime
  retval$mxVersion <- regFit@output$mxVersion
  retval$modelName <- model@name
  plan <- regFit@runstate$compute
  if (is(plan, "MxComputeSequence")) {
    gd <- plan$steps[['GD']]
    if (is(gd, "MxComputeGradientDescent")) {
      retval$optimizerEngine <- gd$engine
    }
  }
  retval$verbose <- verbose
  class(retval) <- "summary.mxmodel"
  return(retval)
}