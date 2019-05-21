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

# This file contains functions that create regularization penalty algebras.

# User-facing Regularize Functions


##' MxRegularizeConstant
##'
##' Constant regularization penalty for no good reason
##'
##' @param name Name of the regularizer object
##' @param penalty How much to penalize
##' 
#' @export 

mxRegularizeConstant <- function(name, penalty) {
  # Note: penalty is not a hyperparameter because we wouldn't want to search over it.
  mxRegularize(c(), how="constant", name=name, otherArgs=list(penalty=penalty))
}

##' MxRegularizeLASSO
##'
##' Least Absolute Selection and Shrinkage Operator regularization
##'
##' @param what A list of parameters to regularize
##' @param name Name of the regularizer object
##' @param lambda strength of the penalty to be applied
##' 
#' @export
mxRegularizeLASSO <- function(what, name, lambda) {
  mxRegularize(what, how="LASSO", hyperparams=list(lambda=lambda), name=name)
}

##' MxRegularizeRidge
##'
##' Ridge regression regularization
##'
##' @param what A list of parameters to regularize
##' @param name Name of the regularizer object
##' @param alpha strength of the penalty to be applied
##' 
#' @export
mxRegularizeRidge <- function(what, name, alpha) {
  mxRegularize(what, how="ridge", hyperparams=list(alpha=alpha), name=name)
}

# Matching algebra-creating functions
# Gotta be a better interface than this.
imxGenerateConstantPenalty <- function(penalty, matrix) {
  return(mxAlgebraFromString(paste(penalty$other$penalty), name=penalty$name))
}

# Just an internal penalty handler for the lasso.
imxGenerateLASSOPenalty <- function(penalty, params, hyperparms) {
  paramName <- params$name
  hyperparmName <- hyperparms$name
  # parameterized to easily allow free parameters or value propagation
  algString <- paste0(hyperparmName, "[1,1] * sum(abs(", paramName, "))")
  return(list(mxAlgebraFromString(algString, name=penalty$name)))
}

# And one for the ridge
imxGenerateRidgePenalty <- function(penalty, matrix, hyperparms) {
  matName <- matrix$name
  # This one's just fixed for now.
  algString <- paste0(hyperparms$name, "[1,1] * t(", matName, ") %*% ", matName)
  return(mxAlgebraFromString(algString, name=penalty$name))
}

# Register regularizers here:
imxRegularizationTypes <- list(LASSO=imxGenerateLASSOPenalty, 
                               ridge=imxGenerateRidgePenalty,
                               constant=imxGenerateConstantPenalty)
