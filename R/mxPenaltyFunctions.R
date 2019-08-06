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
##' @param lambda strength of the penalty to be applied at starting values (default 0)
##' @param lambda.step step function for lambda step (default .01)
##' @param lambda.max end of lambda range (default .4)
##' @param lambda.min minimum lambda value (default lambda)
##' 
#' @export
mxRegularizeLASSO <- function(what, name, lambda=0, lambda.step=.01, lambda.max=NA, lambda.min=NA) {
  if(is.na(lambda.min)) lambda.min=lambda
  if(is.na(lambda.max)) lambda.max=lambda.min+40*lambda.step
  mxRegularize(what, how="LASSO", hyperparams=list(lambda=lambda), hpranges=list(lambda=seq(lambda, lambda.max, by=lambda.step)), name=name)
}

##' MxRegularizeRidge
##'
##' Ridge regression regularization
##'
##' @param what A list of parameters to regularize
##' @param name Name of the regularizer object
##' @param lambda strength of the penalty to be applied at start (default 0)
##' @param lambda.step lambda step during penalty search (default 0.01)
##' @param lambda.max when to end the lambda search (default 0.4) 
##' @param lambda.min minimum lambda value (default lambda)
##' 
#' @export
mxRegularizeRidge <- function(what, name, lambda=0, lambda.step=.01, lambda.max=.4, lambda.min=NA) {
  if(is.na(lambda.min)) lambda.min=lambda
  mxRegularize(what, how="ridge", hyperparams=list(lambda=lambda),
               hprange = list(lambda=seq(lambda, lambda.max, lambda.step)), name=name)
}

##' MxRegularizeElasticNet
##'
##' Elastic net regularization
##'
##' @param what A list of parameters to regularize
##' @param name Name of the regularizer object
##' @param alpha strength of the mixing parameter to be applied at start (default 0.5).  Note that 0 indicates a ridge regression with penalty \deqn{\frac{lambda}{2}}{lambda / 2}, and 1 indicates a LASSO regression with penalty lambda.
##' @param alpha.step alpha step during penalty search (default 0.1)
##' @param alpha.max when to end the alpha search (default 1) 
##' @param lambda strength of the penalty to be applied at starting values (default 0)
##' @param lambda.step step function for lambda step (default .01)
##' @param lambda.max end of lambda range (default .4)
##' @param lambda.min beginning of the lambda range (default lambda)
##' @param alpha.min beginning of the alpha range (default 0)
##' 
##' @details Applies elastic net regularization.  Elastic net is a weighted combination of ridge and LASSO penalties.
##' 
#' @export
mxRegularizeElasticNet <- function(what, name, 
                                   alpha=0,  alpha.step=.1,  alpha.max=1,
                                   lambda=0, lambda.step=.1, lambda.max=.4,
                                   alpha.min=NA, lambda.min=NA) {
  if(is.na(lambda.min)) lambda.min=lambda
  if(is.na(alpha.min)) alpha.min=alpha
  mxRegularize(what, how="elasticNet", hyperparams=list(alpha=alpha, lambda=lambda),
               hprange = list(alpha=seq(alpha, alpha.max, alpha.step), lambda=seq(lambda, lambda.max, by=lambda.step)), name=name)
}

# Matching algebra-creating functions
# Gotta be a better interface than this.
imxGeneratePenaltyConstant <- function(penalty, matrix) {
  return(mxAlgebraFromString(paste(penalty$other$penalty), name=penalty$name))
}

# Just an internal penalty handler for the lasso.
imxGeneratePenaltyLASSO <- function(penalty, params, hyperparms) {
  paramName <- params$name
  hyperparmName <- hyperparms$name
  # parameterized to easily allow free parameters or value propagation
  algString <- paste0(hyperparmName, "[1,1] * sum(abs(", paramName, "))")
  return(list(mxAlgebraFromString(algString, name=penalty$name)))
}

# And one for the ridge
imxGeneratePenaltyRidge <- function(penalty, matrix, hyperparms) {
  matName <- matrix$name
  # This one's just fixed for now.
  algString <- paste0(hyperparms$name, "[1,1] * t(", matName, ") %*% ", matName)
  return(mxAlgebraFromString(algString, name=penalty$name))
}

# And one for the ridge
imxGeneratePenaltyElasticNet <- function(penalty, matrix, hyperparms) {
  matName <- matrix$name
  hyperparmName <- hyperparms$name
  # Lambda * [
  #           (1-alpha) * sqrt(t(vec) %*% vec) +  # ridge
  #           alpha * sum(abs(vec))               # LASSO
  #]
  
  algString <- paste0(hyperparmName, "[1,1] * (",
    "(1-", hyperparmName, "[2,1]) * t(", matName, ") %*% ", matName, ") + ",
    hyperparmName, "[2,1] * sum(abs(", matName, "))" )
  return(mxAlgebraFromString(algString, name=penalty$name))
}

imxGeneratePenaltyAdaptiveLASSO <- function(penalty, matrix, hyperparms) {
  matName <- matrix$name
  # Needs a compute plan to complete.
  algString <- paste0(hyperparms$name, "[1,1] * t(", matName, ") %*% ", matName)
  return(mxAlgebraFromString(algString, name=penalty$name))
}

# Register regularizers here:
imxRegularizationTypes <- list(LASSO=imxGeneratePenaltyLASSO, 
                               ridge=imxGeneratePenaltyRidge,
                               constant=imxGeneratePenaltyConstant,
                               elasticNet=imxGeneratePenaltyElasticNet
                               )

