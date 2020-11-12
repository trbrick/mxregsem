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

# This file contains plotting-focused helper functions.

nameTransform <- function(nameList) {
  nameList <- unlist(nameList)
  return(gsub("[[:punct:]]", ".", nameList))
}

#' Plotting function to see some of the regularization outputs cleanly
#'
#' @param model an MxRegularizedModel object to regularize
#' 
#' @details Plots of the regularization tradeoffs
#' 
#' @return None
#' @import OpenMx
#' @importFrom tidyr gather 
#' @importFrom dplyr %>% select
#' @importFrom ggplot2 ggplot geom_line facet_grid aes_string facet_grid
#' 
#' @export
#' 
plotReg <- function(model) {
  if(length(model$output) == 0 || is.null(model$output$regularizationSearchOutcomes)) {
    stop("plotReg only works on fitted regsem models.")
  }
  penalties <- model@regularizations
  regSearchResult <- model$output$regularizationSearchOutcomes
  for(penalty in penalties) {
    params <- penalty@params
    hparams <- paste(penalty@name, names(penalty@hyperparameters), sep="_")
    if(length(hparams) == 1) {
      meltResult <- regSearchResult %>% 
        dplyr::select(params, hparams, "EBIC") %>%
        tidyr::gather("Variable", "Value", -(hparams))
        print(ggplot2::ggplot(meltResult, ggplot2::aes_string(x = hparams[1], y="Value")) + 
                ggplot2::geom_line() + 
                ggplot2::facet_grid(Variable ~ ., scales = "free_y"))
    } else {
      meltResult <- regSearchResult %>% 
        dplyr::select(params, hparams, "EBIC", "DF", "ll") %>%
        tidyr::gather("Variable", "Value", -(hparams)) 
        for(hp in hparams) {
          print(ggplot2::ggplot(meltResult, ggplot2::aes_string(x = hp, y="Value", color="Variable")) + 
                  ggplot2::geom_line( ggplot2::aes_string(group=setdiff(hparams, hp), color=setdiff(hparams, hp))) + 
                  # ggplot2::geom_hline(ggplot2::aes(xin=)) +
                  ggplot2::facet_grid(Variable ~ ., scales = "free_y"))
                  
        }
      }
  }
  
}
