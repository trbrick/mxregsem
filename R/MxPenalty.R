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

# ##' MxPenalty
# ##'
# ##' This is an internal class and should not be used directly.
# ##'
# ##' @aliases
# ##' $<-,MxPenalty-method
# ##' [[<-,MxPenalty-method
#' @import OpenMx
setClass(Class = "MxPenalty",
          representation = representation(
           name = "character",
           type = "character",
           params = "character",
           hyperparameters = "list",
           other = "list"
          )
         )

setMethod("initialize", "MxPenalty",
          function(.Object, name, type, reg_params, hyperparams=NULL, otherArgs=NULL) {
            .Object@name <- name
            .Object@type <- type
            .Object@params <- reg_params
            .Object@hyperparameters <- hyperparams
            .Object@other <- otherArgs
            return(.Object)
          }
)

setMethod("$", "MxPenalty", imxExtractSlot)

setReplaceMethod("$", "MxPenalty",
                 function(x, name, value) {
                   stop("Error NYI: Can't change penalty values directly.")
                 }
)

setMethod("names", "MxPenalty", slotNames)

##' This function creates a regularization penalty object
##'
##' @param what which parameters or model elements to regularize
##' @param how what kind of regularization function to use.  Currently supported: "lasso", "ridge"
##' @param hyperparams a named list of the hyperparameter starting values. Used for initial setup.
##' @param name the name of the regularization object to be created
##' @param other  Other arguments to the associated regularization functions (see mxRegularizeLASSO and mxRegularizeRidge)
mxRegularize <- function(what, how=names(imxRegularizationTypes), hyperparams=list(), name=NULL, otherArgs=list()) {
    how <- match.arg(how)
    if(is.null(name)) {
      name <- imxUntitledName()
    }
    if(name == "penalty_algebra") {
      stop("NYI: The name 'penalty_algebra' is reserved in regularization.  Sorry, buck-o.")
    }
    new("MxPenalty", name=name, type=how, reg_params=what, hyperparams=hyperparams, otherArgs=otherArgs)
}
