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

# #' Add references levels to beta names created by \code{model.matrix}
# #' 
# #' @return a character vector of refined beta names
# #' @param design a formula specifying the design
# #' @param data a data.frame used for building the model matrix
# #'
# #' @author Yoann Pradat
# #'
# #' @references internal
# refine_beta_names <- function(design, data){
#   data_names <- colnames(data)
#   model_matrix <- model.matrix(design, data)
#   beta_names <- colnames(model_matrix)
# 
#   for (data_name in data_names){
#     data_col <- data[[data_name]] 
#     if (class(data_col) == "factor"){
#       ref_level <- levels(data_col)[1]
#       for (oth_level in levels(data_col)){
#         if (!oth_level==ref_level){
#           old_name <- paste0(data_name, oth_level)
#           new_name <- paste(data_name, oth_level, "vs", ref_level, sep="_")
# 
#           # change name of single var and last term interactions
#           # takes care of ambiguous cases as beta names that are susbtring of other beta names
#           # for instance
#           # - gsub("genotypeII$", "genotype_II_vs_I", x="genotypeIII") will not change x
#           # - gsub("genotypeII$", "genotype_II_vs_I", x="genotypeII") will not change x
#           beta_names <- gsub(paste0(old_name, "$"), new_name, beta_names)
# 
#           # first interactions
#           beta_names <- gsub(paste0("^", old_name, ":"), paste0(new_name, ":"), beta_names)
# 
#           # middle interactions
#           beta_names <- gsub(paste0(":", old_name, ":"), paste0(":", new_name, ":"), beta_names)
#         }
#       }
#     }
#   }
# 
#   beta_names <- gsub(":", ".", beta_names)
#   beta_names <- gsub("\\(|\\)", "", beta_names)
# 
#   beta_names
# }

# refine_beta_names <- function(design, data, main_beta_names=NULL, interacting_beta_names=NULL){
#   design_char <- paste(design)
#   design_char <- strsplit(design_char, "~")[[2]]
#   split_design_names <- unique(trimws(unlist(strsplit(design_char, "\\+|\\-"))))
#   which_inter <- grepl(":|\\*", split_design_names)
#   main_design_names <- unique(split_design_names[!which_inter])
#   inter_design_names <- lapply(split_design_names[which_inter], function(x) trimws(unlist(strsplit(x, "\\*|:"))))
#   inter_design_names <- unique(unlist(inter_design_names))
# 
#   design_names <- unique(main_design_names, inter_design_names)
#   model_matrix <- model.matrix(design, data)
#   beta_names <- colnames(model_matrix)
# 
#   for (design_name in design_names){
#     data_col <- data[[design_name]] 
#     if (class(data_col) == "factor"){
#       ref_level <- levels(data_col)[1]
#       for (oth_level in levels(data_col)){
#         if (!oth_level==ref_level){
#           # for genotype + genotype:condition design, do not update genotype names in the interaction part
#           # as there will be one interaction term for each level of genotype
# 
#           # for condition + genotype:condition design, do not update condition names in the interaction part
#           # as there will be one interaction term for each level of condition 
# 
#           # for condition + genotype + genotype:condition design, update both condition and genotype names in the 
#           # interaction part
# 
#           old_name <- paste0(design_name, oth_level)
#           new_name <- paste(design_name, oth_level, "vs", ref_level, sep="_")
# 
#           # change name of single var and last term interactions
#           # takes care of ambiguous cases as beta names that are susbtring of other beta names
#           # for instance
#           # - gsub("genotypeII$", "genotype_II_vs_I", x="genotypeIII") will not change x
#           # - gsub("genotypeII$", "genotype_II_vs_I", x="genotypeII") will not change x
#           beta_names <- gsub(paste0(old_name, "$"), new_name, beta_names)
# 
#           # first interactions
#           beta_names <- gsub(paste0("^", old_name, ":"), paste0(new_name, ":"), beta_names)
# 
#           # middle interactions
#           beta_names <- gsub(paste0(":", old_name, ":"), paste0(":", new_name, ":"), beta_names)
#           }
#         }
#       }
#     }
#   }
# 
#   beta_names <- gsub(":", ".", beta_names)
#   beta_names <- gsub("\\(|\\)", "", beta_names)
# 
#   beta_names
# }


#' Produce a contrast vector
#' 
#' @return a character vector of refined beta names
#' @param contrast a character vector specifying the contrast beween refined beta names (see
#' \code{\link{refine_beta_names}}
#' @param design a formula specifying the design
#' @param data a data.frame used for building the model matrix
#'
#' @importFrom limma makeContrasts
#' @importFrom stats model.matrix
#'
#' @author Yoann Pradat
#'
#' @references internal
get_contrast_vector <- function(contrast, design, data){
  #beta_names <- refine_beta_names(design,data)

  data_names <- colnames(data)
  model_matrix <- model.matrix(design, data)
  beta_names <- colnames(model_matrix)
  beta_names <- gsub(":", ".", beta_names)
  beta_names <- gsub("\\(|\\)", "", beta_names)

  cmd <- paste("contrast_vec <- makeContrasts(", contrast, ", levels = beta_names)", sep = '"')
  eval(parse(text = cmd))
  contrast_vec
}


#' Build a list of character contrasts of B-A for an interacting variable.
#' 
#' @param levels_cond the factor levels of the condition variable
#' @param levels_inte the factor levels of the interaction variable
#' @param cond_name the name of the condition variable
#' @param inte_name the name of the interaction variable
#' @param level_B_cond the B level to be contrasted with A
#' @param level_A_cond the A level to be contrasted with B
#' @param cond_main_effect is the condition a main effect
#' @param inte_main_effect is the interaction a main effect
#'
#' @author Yoann Pradat
#'
#' @export
get_contrast_condition_B_minus_A_factors_interaction <- function(levels_cond, levels_inte, cond_name="condition",
                                                                 inte_name="interaction",
                                                                 level_B_cond=NULL, 
                                                                 level_A_cond=NULL, 
                                                                 cond_main_effect=T, inte_main_effect=F){

  if (!cond_main_effect & !inte_main_effect){
    stop("when using an interaction between two factor variables, at least one variable should be a main effect")
  }

  ref_cond <- levels_cond[1]
  ref_inte <- levels_inte[1]

  if (is.null(level_B_cond)){
    message("unspecified level condition B, choosing the last level")
    level_B_cond <- levels_cond[length(levels_cond)]
  }

  if (is.null(level_A_cond)){
    message("unspecified level condition A, choosing the first (i.e the reference) level")
    level_A_cond <- ref_cond
  }

  if (level_A_cond == level_B_cond){
    stop("condition levels A and B to be contrasted are equal")
  }

  contrasts <- list()

  if (cond_main_effect & !inte_main_effect){
    for (level_inte in levels_inte){
      contrast_name <- paste(cond_name, level_B_cond, "vs", level_A_cond, inte_name, level_inte, sep="_")

      if (level_inte == ref_inte){
        if (level_A_cond == ref_cond){
          contrast_expr <- paste0(cond_name, level_B_cond)
        } else if (level_B_cond == ref_cond) {
          contrast_expr <- paste0("-", cond_name, level_A_cond)
        } else {
          contrast_expr <- paste0(cond_name, level_B_cond, "-", cond_name, level_A_cond)
        }
      } else {
        contrast_expr <- paste(paste0(cond_name, level_B_cond, ".", inte_name, level_inte),
                               paste0(cond_name, level_A_cond, ".", inte_name, level_inte),
                               sep="-")
      }
      contrasts[[contrast_name]] <- contrast_expr
    }
  } else {
    stop ("Not yet implemented")
  }

  contrasts
}
