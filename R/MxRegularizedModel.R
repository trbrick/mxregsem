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

#' @import OpenMx
# #' MxRegularizedModel
# ##'
# ##' This is an internal class and should not be used directly.
# ##'
# ##' @aliases
# ##' $<-,MxRegularizedModel-method
# ##' [[<-,MxRegularizedModel-method
setClass(Class = "MxRegularizedModel",
         slots = representation(
           regularizations = "list"
         ),
         contains = "MxModel")

imxInitRegularizedModel <- function(submodel1) {
            model <- new("MxRegularizedModel")
            model@submodels[[submodel1$name]] <- submodel1
            if(is.null(model@algebras[['penaltyFunction']])) {
              model@algebras[['penaltyFunction']] <- mxAlgebraFromString(paste0(submodel1$name,".fitfunction"), name="penaltyFunction")
            }
            if(is.null(model@fitfunction)) {
              model@fitfunction <- mxFitFunctionAlgebra("penaltyFunction")
            }
            model@name <- paste0(submodel1$name, "_Reg")
            return(model)
          }

setMethod("imxVerifyModel", "MxRegularizedModel",
          function(model) {
            if (length(model@submodels) == 0) {
              msg <- paste("The regularized model", omxQuotes(model@name),
                           "does not regularize anything.")
              stop(msg, call. = FALSE)
            }
            fitfunction <- model@fitfunction
            if (is.null(fitfunction) || !is(fitfunction, "MxFitFunctionAlgebra")) {
                  msg <- paste("The regularized model", omxQuotes(model@name),
                               "is somehow malformed. Sorry!")
                  stop(msg, call. = FALSE)
            }
            if (length(model@submodels) > 0) {
              return(all(sapply(model@submodels, imxVerifyModel)))
            }
            return(TRUE)
          })

setMethod("print", "MxRegularizedModel", function(x,...) { 
  displayRegModel(x, TRUE) 
})



#' Create a (regularizable) mxModel
#'
#' Identical to OpenMx::mxModel, except the result is an MxRegularizedModel.
#'
#' @importFrom stats sd
#' @export
#' 
mxModel <- function(model = NA, ..., manifestVars = NA, latentVars = NA,
                              remove = FALSE, independent = NA, type = NA, name = NA) {

  retval <- regFirstArgument(model, name)  # Reimplemented `firstArgument()` function from MxModel
  first <- retval[[1]]
  regModel <- retval[[2]]
  name  <- retval[[3]]
  regModel@submodels[[1]] <- OpenMx:::typeArgument(regModel@submodels[[1]], type)
  if(remove==TRUE) {
    # Removal punts upstream.  Use with care.
    return(OpenMx::mxModel(regModel, first, ..., name=name, manifestVars=manifestVars,
                          latentVars=latentVars, remove=remove, independent=independent, type=type))
  }
  lst <- c(first, list(...))
  lst <- unlist(lst)
  filter <- sapply(lst, is, "MxModel")
  submodels <- lst[filter]
  lst <- lst[!filter]
  filter <- sapply(lst, is, "MxPenalty")
  if(any(filter)) {
    regModel <- mxRegularizedAddPenalty(regModel, lst[filter])
  }
  lst <- lst[!filter]
  regModel@submodels[[1]] <- imxModelBuilder(regModel@submodels[[1]], lst, name, manifestVars,
                           latentVars, submodels, remove, independent)
  # Handle renaming
  regModel <- regUpdateNames(regModel)
  
  return(regModel)
}

regUpdateNames <- function(regModel) {
  if(!is(regModel, "MxRegularizedModel")) {
    return(regModel)
  }
  
  # Handle renaming
  regModel@name <- paste0(regModel@submodels[[1]]@name, "_Reg")
  
  # Namespace maintenance.  Note: use mxRename to really rename.
  names(regModel@submodels)[[1]] <- regModel@submodels[[1]]@name
  
  return(regModel)
}

# Helpers to help build regularized model objects
regFirstArgument <- function(model, name) {
  first <- NULL
  defaultType <- imxModelTypes[[getOption("mxDefaultType")]]
  
  # Check if it's a regularized model already:
  if(is(model, "MxRegularizedModel")) {
    regModel <- model
  } else if (is(model, "MxModel")) {
    regModel <- imxInitRegularizedModel(model)
  } else {
    # Creating a new model: need to make both model and regmodel here.
    if (OpenMx:::single.na(model)) {
      # Nothing doing here.
    } else if (typeof(model) == "character") {
      name <- model
    } else if (isS4(model)) {
      first <- model
    } else {
      first <- list(model)
    }
    if (length(name) > 0 && is.na(name)) {
      name <- imxUntitledName()
    }
    imxVerifyName(name, -1)
    model <- new(defaultType, name)
    regModel <- imxInitRegularizedModel(model)
  }
  return(list(first, regModel, name))
}

mxRegularizedAddPenalty <- function(regModel, penalties) {
  pNames <- sapply(penalties, getElement, "name")
  regNames <- sapply(regModel@regularizations, getElement, "name")
  conflicts <- pNames %in% regNames
  if(any(conflicts)) {
    for(aMatch in which(conflicts)) {
      # Remove original regularizer; replace below.
      regModel@regularizations[[aMatch]] <- NULL
      if(!is.null(regModel@matrices[[paste0(aMatch, "_Params")]])) {
        regModel@matrices[[paste0(aMatch, "_Params")]] <- NULL
      }
      if(!is.null(regModel@matrices[[paste0(aMatch, "_Penalty")]])) {
        regModel@matrices[[paste0(aMatch, "_Penalty")]] <- NULL
      }
      regModel@algebras[[aMatch]] <- NULL
    }
  }
  
  paramVals <- omxGetParameters(regModel, fetch="values")
  # Add algebras for each new penalty:
  for(aPenalty in penalties) {
    # Register penalty
    regModel@regularizations[[aPenalty$name]] <- aPenalty
    
    # Generate parameter vectors and create the penalty matrix
    reg_params <- aPenalty$params
    penalty_Mat <- NA
    if(!is.null(reg_params) && length(reg_params) > 0) {
      pVals <- paramVals[reg_params]
      paramsFree <- rep(TRUE, length(pVals))
      paramsFree[grep("[\\[]", reg_params)] <- FALSE
    
      # Doesn't handle def vars yet.
      if(length(grep("data[\\.]", reg_params))>0) {
        stop("Error NYI in regularization penalty ", aPenalty$name,
              ": cannot regularize definition vars.")
      }
      
      # Add penalty to the regModel
      paramsMat <- mxMatrix("Full", nrow=length(paramsFree), 
                             ncol=1, values=pVals, 
                             free=paramsFree, labels=reg_params, 
                             name=paste0(aPenalty$name, "_Params"))
      regModel <- OpenMx::mxModel(regModel, paramsMat)
      # Add hyperparameters
      if(length(aPenalty@hyperparameters) > 0) {
        penaltyMat <- mxMatrix("Full", nrow=length(aPenalty@hyperparameters), 
                             ncol=1, free=FALSE, 
                             values=unlist(aPenalty@hyperparameters),
                             name=paste0(aPenalty$name, "_Penalty"),
                             labels=paste(aPenalty@name, 
                                          names(aPenalty@hyperparameters), 
                                          sep="_"),
                            )
        regModel <- OpenMx::mxModel(regModel, penaltyMat)
      }
    }
    
    regModel <- OpenMx::mxModel(regModel, 
        imxRegularizationTypes[[aPenalty@type]](aPenalty, paramsMat, penaltyMat))
  }
  
  # Rebuild global penalty algebra:
  penaltyString <- paste0(regModel@submodels[[1]]$name, ".fitfunction")
  if(!is.null(regModel@regularizations) && length(regModel@regularizations) != 0) {
    penaltyString <- paste(penaltyString, names(regModel@regularizations), sep="+")
  }
  
  regModel@algebras[["penaltyFunction"]] <- mxAlgebraFromString(penaltyString, name="penaltyFunction")
  
  return(regModel)
}

#' Display a regularized model:
displayRegModel <- function(regmodel, expand = FALSE) {
  omxQuotes <- OpenMx::omxQuotes # Better done with namespace package.
  model <- regmodel@submodels[[1]]
  cat("MxModel", omxQuotes(model@name), '(Regularized)\n')
  cat("type :", OpenMx::imxTypeName(model), '\n')
  cat("$matrices :", omxQuotes(names(model@matrices)), '\n')
  cat("$algebras :", omxQuotes(names(model@algebras)), '\n')
  cat("$constraints :", omxQuotes(names(model@constraints)), '\n')
  cat("$intervals :", omxQuotes(names(model@intervals)), '\n')
  cat("$regularizers :", omxQuotes(names(regmodel@regularizers)), '\n')
  
  # latentVars and manifestVars should really be considered
  # an implementation detail of RAM and LISREL type models with paths.
  # We currently do not return anything when an attempt
  # is made to access these slots using $ notation.
  # The proper thing to do is probably to ignore the @ slots
  # and extract the variable information from the RAM F matrix dimnames
  # or equivalent LISREL matrix because the slots are only
  # used by mxPath. Matrix constructed models do not use the @ slots.
  
  if (length(model@latentVars) == 0 || OpenMx::imxTypeName(model) %in% "default") {
    cat("$latentVars : none\n")
  } else if (is.character(model@latentVars)) {
    cat("$latentVars :", omxQuotes(model@latentVars), '\n')
  } else {
    cat("$latentVars :\n")
    print(format(model@latentVars))
  }
  if (length(model@manifestVars) == 0 || imxTypeName(model) %in% "default") {
    cat("$manifestVars : none\n")
  } else if (is.character(model@manifestVars)) {
    cat("$manifestVars :", omxQuotes(model@manifestVars), '\n')
  } else {
    cat("$manifestVars :\n")
    print(format(model@manifestVars))
  }
  data <- model@data
  if (is.null(data)) {
    cat("$data : NULL\n")
  } else {
    if (is(data, "MxDataDynamic")) {
      cat("$data type:", omxQuotes(data@type), '\n')
      cat("$data$expectation :", omxQuotes(data@expectation), "\n")
    } else {
      cat("$data :", nrow(data@observed), 
          "x", ncol(data@observed), "\n")
      if(length(data@means) == 1 && is.na(data@means)) {
        cat("$data means : NA\n")
      } else {
        cat("$data means : 1 x", length(data@means), "\n")
      }
      cat("$data type:", omxQuotes(data@type), '\n')
    }
  }
  cat("$submodels :", omxQuotes(names(model@submodels)), '\n')
  expectation <- model@expectation
  fitfunction <- model@fitfunction
  regularization <- regmodel@regularization
  compute <- model@compute
  if (is.null(expectation)) {
    expectationType <- "NULL"
  } else {
    expectationType <- class(expectation)[[1]]
  }
  cat("$expectation :", expectationType, '\n')
  
  if (is.null(fitfunction)) {
    fitfunctionType <- "NULL"
  } else {
    fitfunctionType <- class(fitfunction)[[1]]
  }
  cat("$fitfunction :", fitfunctionType, '\n')
  
  if (is.null(compute)) {
    computeType <- "NULL"
  } else {
    computeType <- class(compute)[[1]]
  }
  cat("$compute :", computeType, '\n')
  cat("$independent :", model@independent, '\n')
  cat("$options :", OpenMx::printOptions(model@options), '\n')
  cat("$output :", length(model@output) > 0, '\n')
  if(is(regmodel, "MxRegularizedModel")) {
    if (is.null(regularization) || length(regularization) == 0) {
      regularizationType <- "NULL"
    } else {
      regularizationType <- paste(sapply(regularization, getElement, "type"), sep=" & ")
    }
    cat("$regularization :", regularizationType)
    if(length(regmodel@matrices) > 0) {
      cat("\n--------REGULARIZATION MATRICES--------\n")
      lapply(regmodel@matrices, print)
    }
    if(length(regmodel@algebras) > 0) {
      cat("\n--------REGULARIZATION ALGEBRAS--------\n")
      lapply(regmodel@algebras, print)
    }
    if(!is.null(regmodel@submodels[[1]])) {
      if (is.null(regmodel@submodels[[1]]@expectation)) {
        regularizationType <- "NULL"
      } else {
        expectationType <- class(regmodel@submodels[[1]]@expectation)[[1]]
      }
      cat("\t-regularized expectation :", expectationType, '\n')
    
      if (is.null(regmodel@submodels[[1]]@fitfunction)) {
        fitfunctionType <- "NULL"
      } else {
        fitfunctionType <- class(regmodel@submodels[[1]]@fitfunction)[[1]]
      }
      cat("\t-regularized fitfunction :", fitfunctionType, '\n')
    }
    
  }
  if(expand) {
    if(length(model@matrices) > 0) {
      cat("\n--------MATRICES--------\n")
      lapply(model@matrices, print)
    }
    if(length(model@algebras) > 0) {
      cat("\n--------ALGEBRAS--------\n")
      lapply(model@algebras, print)
    }
    if(length(model@constraints) > 0) {
      cat("\n--------CONSTRAINTS--------\n")
      lapply(model@constraints, print)
    }
    if(!is.null(model@data) > 0) {
      cat("\n--------DATA--------\n")
      print(model@data)
    }
    if(!is.null(model@expectation) > 0) {
      cat("\n--------EXPECTATION FUNCTION--------\n")
      print(model@expectation)
    }
    if(!is.null(model@fitfunction) > 0) {
      cat("\n--------FIT FUNCTION--------\n")
      print(model@fitfunction)
    }
    if(!is.null(model@compute) > 0) {
      cat("\n--------COMPUTE--------\n")
      print(model@compute)
    }
    if(length(model@output) > 0) {
      cat("\n--------OUTPUT--------\n")
      print(model@output)
    }
    if(length(model@submodels) > 0) {
      cat("\n--------SUBMODELS--------\n")
      lapply(model@submodels, print)
    }
    if(length(model@options) > 0) {
      cat("\n--------OPTIONS--------\n")
      print(model@options)
    }
  }
  invisible(model)
}

# Overwrite renaming:
##' Not documented.  See OpenMx::mxRename.
##' @export
mxRename <- function(model, newname, oldname = NA) {
  model <- OpenMx::mxRename(model, newname, oldname)
  # Apply namespace corrections.
  model <- regUpdateNames(model)
  return(model)
}

parseLikelihoodArg <- function(input, arg) {
  input <- input[[arg]]
  if (is.null(input)) {
    return(input)
  } else if (is.numeric(input)) {
    return(input)
  } else if (is(input, "MxModel")) {
    name <- input@name
    if (is.null(input@fitfunction)) {
      stop(paste(omxQuotes(name), "model passed",
                 "to summary function does not",
                 "have top-level fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    if (length(input@fitfunction@result) != 1) {
      stop(paste(omxQuotes(name), "model passed to summary",
                 "function does not have a 1x1 matrix",
                 "result in fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    return(input@fitfunction@result[1,1])
  } else if(is.list(input) && length(input)==2) {
    stop(paste("List of length two (illegal argument) passed to", omxQuotes(arg),
               "argument of summary function. You probably meant to use",
               "the refModels argument instead."), call. = FALSE)
  } else {
    stop(paste("Illegal argument passed to", omxQuotes(arg),
               "argument of summary function in",
               deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
  }
}

parseDfArg <- function(input, arg) {
  input <- input[[arg]]
  if (is.null(input)) {
    return(input)
  } else if (is.numeric(input)) {
    return(input)
  } else if (is(input, "MxModel")) {
    name <- input@name
    if (is.null(input@fitfunction)) {
      stop(paste(omxQuotes(name), "model passed",
                 "to summary function does not",
                 "have top-level fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    if (length(input@fitfunction@result) != 1) {
      stop(paste(omxQuotes(name), "model passed to summary",
                 "function does not have a 1x1 matrix",
                 "result in fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    return(summary(input)$degreesOfFreedom)
  } else {
    stop(paste("Illegal argument passed to", omxQuotes(arg),
               "argument of summary function in",
               deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
  }
}

parseLikelihoodArg <- function(input, arg) {
  input <- input[[arg]]
  if (is.null(input)) {
    return(input)
  } else if (is.numeric(input)) {
    return(input)
  } else if (is(input, "MxModel")) {
    name <- input@name
    if (is.null(input@fitfunction)) {
      stop(paste(omxQuotes(name), "model passed",
                 "to summary function does not",
                 "have top-level fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    if (length(input@fitfunction@result) != 1) {
      stop(paste(omxQuotes(name), "model passed to summary",
                 "function does not have a 1x1 matrix",
                 "result in fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    return(input@fitfunction@result[1,1])
  } else if(is.list(input) && length(input)==2) {
    stop(paste("List of length two (illegal argument) passed to", omxQuotes(arg),
               "argument of summary function. You probably meant to use",
               "the refModels argument instead."), call. = FALSE)
  } else {
    stop(paste("Illegal argument passed to", omxQuotes(arg),
               "argument of summary function in",
               deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
  }
}

parseDfArg <- function(input, arg) {
  input <- input[[arg]]
  if (is.null(input)) {
    return(input)
  } else if (is.numeric(input)) {
    return(input)
  } else if (is(input, "MxModel")) {
    name <- input@name
    if (is.null(input@fitfunction)) {
      stop(paste(omxQuotes(name), "model passed",
                 "to summary function does not",
                 "have top-level fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    if (length(input@fitfunction@result) != 1) {
      stop(paste(omxQuotes(name), "model passed to summary",
                 "function does not have a 1x1 matrix",
                 "result in fitfunction in",
                 deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
    }
    return(summary(input)$degreesOfFreedom)
  } else {
    stop(paste("Illegal argument passed to", omxQuotes(arg),
               "argument of summary function in",
               deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
  }
}

refToLikelihood <- function(model) {
  if (is(model, "MxModel")) {
    if (!model@.wasRun) stop("Reference model must be run to obtain fit indices")
    model$output$Minus2LogLikelihood
  } else if (is.list(model)) {
    model[[1]]
  } else {
    stop(paste("Illegal argument passed to refModels",
               "argument of summary function in",
               deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
  }
}

refToDof <- function(model) {
  if (is(model, "MxModel")) {
    if (!model@.wasRun) stop("Reference model must be run to obtain fit indices")
    return(summary(model)$degreesOfFreedom)
  } else if (is.list(model)) {
    model[[2]]
  } else {
    stop(paste("Illegal argument passed to refModels",
               "argument of summary function in",
               deparse(width.cutoff = 400L, sys.call(-1))), call. = FALSE)
  }
}

#' computeOptimizationStatistics
#'
#' This is overridden from OpenMx.  This function is exported for people
#'    who already know what they are doing.
#'
computeOptimizationStatistics <- function(model, numStats, useSubmodels,
                                          saturatedDoF, independenceDoF,
                                          retval) {
  # get estimated parameters
  estimates <- model@output$estimate
  # should saturated/independence models include means?
  if(length(model@runstate$datalist)==1){
    type <- model@runstate$datalist[[1]]@type
    means <- model@runstate$datalist[[1]]@means
    # if there's raw data, then use means in saturated/independence models
    if(type=="raw"){
      useMeans <- TRUE
    } else {
      # if there's not raw data, only use means if they're present
      if((dim(means)[2]==1)&is.na(means[1,1])){
        useMeans <- FALSE
      } else{
        useMeans <- TRUE
      }
    }
    # number of variables
    if(model@runstate$datalist[[1]]@type != 'raw'){
      nvar <- dim(model@runstate$datalist[[1]]@observed)[2]
    } else if( length(model@runstate$expectations) == 1 ) {
      nvar <- length(model@runstate$expectations[[1]]@dims)
    } else {
      nvar <- 0
    }
    # if there are multiple or zero datalists, then do nothing
  } else {
    useMeans <- NA
    nvar <- 0
  }
  # how many thresholds does each variable have (needed for saturated and independence DoF calculation)
  # grab the expectation
  obj <- model@runstate$expectation
  # grab the thresholdLevels object and expected means; punt if there is more than one expectation
  if (length(obj)==1){
    if ("thresholdLevels" %in% slotNames(obj[[1]])){
      thresholdLevels <- obj[[1]]@thresholdLevels
      if (length(thresholdLevels)==0){thresholdLevels <- rep(NA, nvar)}
    } else {
      thresholdLevels <- rep(NA, nvar)
    }
  } else {
    thresholdLevels <- NULL
  }
  # number of continuous variables, provided there is just one expectation
  if (!is.null(thresholdLevels)){
    continuous <- sum(is.na(thresholdLevels))
  } else{
    continuous <- NA
  }
  # number of thresholds in the model
  if (!is.null(thresholdLevels)){
    thresh <- sum(thresholdLevels, na.rm=TRUE)
  } else{
    thresh <- NA
  }
  # constraints, parameters, model degrees of freedom
  retval[['constraints']] <- OpenMx:::calculateConstraints(model, useSubmodels)
  retval[['estimatedParameters']] <- nrow(retval$parameters)
  if(any(sapply(obj,function(x){"MxExpectationGREML" %in% class(x)}))){
    retval[['estimatedParameters']] <- retval[['estimatedParameters']] +
      sum(sapply(obj,OpenMx:::imxExtractSlot,name="numFixEff"))
  }
  if(!is.null(retval[['regularizedParameters']])) {
    retval$estimatedParameters <- retval$estimatedParameters-retval$regularizedParameters
  }
  # browser()
  # TODO: Intervene here!
  if (is.null(numStats)) {
    retval[['observedStatistics']] <- OpenMx:::observedStatistics(model, useSubmodels, sum(retval$constraints))
  } else {
    retval[['observedStatistics']] <- numStats
  }
  retval[['degreesOfFreedom']] <- retval$observedStatistics - retval$estimatedParameters
  # calculate or populate saturated degrees of freedom
  if(is.null(saturatedDoF)) {
    retval[['saturatedDoF']] <- retval$observedStatistics - (nvar * (nvar-1) / 2 + continuous*(1+useMeans) + thresh)
  } else {
    retval[['saturatedDoF']] <- saturatedDoF
  }
  #The "saturated model" has no sensible definiton with GREML expectation:
  if(any(sapply(obj,function(x){"MxExpectationGREML" %in% class(x)}))){
    retval[['saturatedDoF']] <- NA
  }
  # calculate or populate independence degrees of freedom
  if(is.null(independenceDoF)) {
    if(!any(sapply(obj,function(x){"MxExpectationGREML" %in% class(x)}))){
      # indDoF = 1 df per continuous variable variance + 1 df per continuous mean + 1 df per threshold
      retval[['independenceDoF']] <- retval$observedStatistics - (continuous*(1+useMeans) + thresh)
    } else{
      #TODO: the GREML expectation doesn't currently have a way to know how many phenotypes there are in every case.
      #For now, leave the GREML independence model undefined
      # #With GREML expectation, the independence model has a variance for each phenotype, and the same fixed effects as the fitted model:
      # retval[['independenceDoF']] <-
      # 	retval$observedStatistics - sum(sapply(obj,function(x){length(x@yvars)})) - sum(sapply(obj,imxExtractSlot,name="numFixEff"))
      retval[['independenceDoF']] <- NA
    }
  } else {
    retval[['independenceDoF']] <- independenceDoF
  }
  # set NULLs to NAs
  if (is.null(retval$saturatedDoF)) {
    retval$SaturatedDoF <- NA
  }
  if (is.null(retval$independenceDoF)) {
    retval$IndependenceDoF <- NA
  }
  retval[['saturatedParameters']] <- retval[['observedStatistics']] - retval[['saturatedDoF']]
  retval[['independenceParameters']] <- retval[['observedStatistics']] - retval[['independenceDoF']]
  # calculate fit statistics
  retval <- OpenMx:::fitStatistics(model, useSubmodels, retval)
  return(retval)
}

#' Summarize an MxRegularizedModel
#'
#' Create a summary of an MxRegularizedModel.
#'
#' @param regFit the regularized model to summarize
#' @param ... values passed to mxSummarize
#' @param verbose gives more summary feedback
#' @param epsilon how close to zero is close enough? Regularization limits for DF correction
#' @return Summary object summarizing the MxRegularizedModel
#' @import OpenMx
#' @export
setMethod("summary", "MxRegularizedModel", function(object, ..., verbose=FALSE, computeRefs=FALSE, epsilon = 1e-6) {
  # summary.MxRegularizedModel <- function(object, ..., verbose=FALSE, epsilon = 1e-6) {
  regFit <- object
  # default penalty_function ("guess") looks for a "lasso" object in the container model.
  dotArguments <- list(...)

  if(length(regFit$submodels) == 1) {
    model <- regFit$submodels[[1]]
  } else {
    model <- regFit
  }
  if (!is.null(dotArguments[["refModels"]])) {
    refModels <- dotArguments[["refModels"]]
    satModel <- refModels[['Saturated']]
    indModel <- refModels[['Independence']]
    saturatedLikelihood <- refToLikelihood(satModel)
    saturatedDoF <- refToDof(satModel)
    independenceLikelihood <- refToLikelihood(indModel)
    independenceDoF <- refToDof(indModel)
  } else {
    saturatedLikelihood <- parseLikelihoodArg(dotArguments, "SaturatedLikelihood")
    saturatedDoF <- parseDfArg(dotArguments, "SaturatedDoF")
    independenceLikelihood <- parseLikelihoodArg(dotArguments, "IndependenceLikelihood")
    independenceDoF <- parseDfArg(dotArguments, "IndependenceDoF")
    if(computeRefs && (is.null(saturatedDoF) || is.null(independenceDoF))) {
      refModels <- mxRefModels(regFit$submodels[[1]], run=FALSE)
      satModel <- OpenMx::mxRun(refModels[['Saturated']], silent=TRUE)
      indModel <- OpenMx::mxRun(refModels[['Independence']], silent=TRUE)
      saturatedLikelihood <- refToLikelihood(satModel)
      saturatedDoF <- refToDof(satModel)
      independenceLikelihood <- refToLikelihood(indModel)
      independenceDoF <- refToDof(indModel)
    }
  }
  
  numObs <- dotArguments$numObs
  numStats <- dotArguments$numStats
  
  # Handle regularization elements
  regList <- c()
  for(penalty in regFit@regularizations) {
    # TODO: Manage uniqueness in a better way
    regged <- omxGetParameters(regFit, free=TRUE, labels=penalty@params)
    zeroedOut <- names(regged)[abs(regged) < epsilon]
    regList <- union(regList, zeroedOut)
  }
  # Win back free parameters for any that's less than zero.
  regOffset <- length(regList)
  
  useSubmodels <- dotArguments$indep
  if (is.null(useSubmodels)) { useSubmodels <- TRUE }
  retval <- list(wasRun=model@.wasRun, stale=model@.modifiedSinceRun)
  retval$parameters <- OpenMx:::parameterList(regFit, useSubmodels)
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
  retval <- OpenMx:::setLikelihoods(regFit, saturatedLikelihood, independenceLikelihood, retval)
  retval[["Minus2LogLikelihood"]] <- mxEvalByName(paste0(model$name, ".fitfunction"), regFit)[1,1]
  retval <- OpenMx:::setNumberObservations(numObs, regFit@runstate$datalist, regFit@runstate$fitfunctions, retval)
  retval[["regularizedParameters"]] <- regOffset
  retval[["zeroedParameters"]] <- regList
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
})
