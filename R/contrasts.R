# @modified: 17 Nov 2020
# @created: 17 Nov 2020
# @author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
# 
#     Institut Gustave Roussy
#     Prism Center
#     114 rue Edouard Vaillant, Villejuif, 94800 France
# 
# Functions for preparing contrasts for each differential analysis procedure.

#' Add references levels to beta names created by \code{model.matrix}
#' 
#' @return a character vector of refined beta names
#' @param design a formula specifying the design
#' @param data a data.frame used for building the model matrix
#'
#' @author Yoann Pradat
#'
#' @references internal
refine_beta_names <- function(design, data){
  data_names <- colnames(data)
  model_matrix <- model.matrix(design, data)
  beta_names <- colnames(model_matrix)

  for (data_name in data_names){
    data_col <- data[[data_name]] 
    if (class(data_col) == "factor"){
      ref_level <- levels(data_col)[1]
      for (oth_level in levels(data_col)){
        if (!oth_level==ref_level){
          old_name <- paste0(data_name, oth_level)
          new_name <- paste(data_name, oth_level, "vs", ref_level, sep="_")

          # change name of single var and last term interactions
          # takes care of ambiguous cases as beta names that are susbtring of other beta names
          # for instance
          # - gsub("genotypeII$", "genotype_II_vs_I", x="genotypeIII") will not change x
          # - gsub("genotypeII$", "genotype_II_vs_I", x="genotypeII") will not change x
          beta_names <- gsub(paste0(old_name, "$"), new_name, beta_names)

          # first interactions
          beta_names <- gsub(paste0("^", old_name, ":"), paste0(new_name, ":"), beta_names)

          # middle interactions
          beta_names <- gsub(paste0(":", old_name, ":"), paste0(":", new_name, ":"), beta_names)
        }
      }
    }
  }

  beta_names <- gsub(":", ".", beta_names)
  beta_names <- gsub("\\(|\\)", "", beta_names)

  beta_names
}

#' Produce a contrast vector
#' 
#' @return a character vector of refined beta names
#' @param contrast a character vector specifying the contrast beween refined beta names (see
#' \code{\link{refine_beta_names}}
#' @param design a formula specifying the design
#' @param data a data.frame used for building the model matrix
#'
#' @importFrom limma makeContrasts
#'
#' @author Yoann Pradat
#'
#' @references internal
get_contrast_vector <- function(contrast, design, data){
  levels <- refine_beta_names(design,data)
  cmd <- paste("contrast_vec <- makeContrasts(", contrast, ", levels = levels)", sep = '"')
  eval(parse(text = cmd))
  contrast_vec
}
