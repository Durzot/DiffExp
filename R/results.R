# @modified: 12 Nov 2020
# @created: 12 Nov 2020
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
# Functions for saving table of results as dataframes with specific update and merge functionnalities.

#' Build table of results
#' 
#' @return a data.frame
#' @param data a named list of columns. The names sould include "padj" if \code{only_significant} is set to TRUE.
#' @param tags a named list of scalar that will appended as constant columns
#' @param only_significant whether only significant results should be recorded
#' @param alpha significance level. Used only if \code{only_signifcant} is set to TRUE.
#' @param df_results_all Preexisting table of results to which the new results will be appended
#' @param df_row Extra row information that is to be appended to the table of results
#' @param file set to NULL to not save onto the disk
#'
#' @author Yoann Pradat
#' @keywords internal
make_table_results <- function(data, tags, only_significant=T, alpha=0.1, df_results_all=NULL, df_row=NULL, 
                               file=NULL){
    # make it a dataframe
    df_results <- data.frame(data)
    extra_tags <- list(time=as.character(Sys.time()))
    df_results <- cbind.data.frame(df_row, c(tags, extra_tags), df_results)

    # decide if all the table is recorded or not
    if (only_significant){
      keep <- !is.na(data$padj) & data$padj < alpha
    } else {
      keep <- rep(T, nrow(df_results))
    }
    df_results <- df_results[keep,]

    # save if specified
    if (!is.null(file)){
      save_update_table(file, df_results, tags)
    }

    # append
    rbind(df_results_all, df_results)
}
 

#' Update an existing table with a new table.
#' 
#' @return a data.frame
#' @param tab a data.frame
#' @param tab_new a data.frame
#' @param tags named list of tags values to decide between merging and updating
#'
#' @author Yoann Pradat
#' @keywords internal
update_table <- function(tab, tab_new, tags){
  if (!all(names(tags) %in% names(tab))){
    stop("specified tags not in existing table")
  }

  mask <- rep(T, nrow(tab))
  for (tag_name in names(tags)){
    mask <- mask & tab[[tag_name]] == tags[tag_name]
  }

  if (sum(mask)>0){
    if (sum(mask)!=nrow(tab_new)){
      stop("the part of table to be updated and the new table do not have the same number of rows")
    } else {
      tab <- tab[!mask,]
    }
  }
  tab <- merge(tab, tab_new, all=T)

  tab
}

#' Update an existing table with a new table.
#' 
#' @return a data.frame
#' @param file path to where the table should be saved 
#' @param tab_new a data.frame
#' @param tags named list of tags values to decide between merging and updating if a table already exists at file
#'
#' @author Yoann Pradat
#' @keywords internal
save_update_table <- function(file, tab_new, tags, verbose=T){
  if (!file.exists(file)){
    if (verbose){
      cat(paste("wrote table at", file), "\n")
    }
  } else {
    tab <- read.table(file, header=T, sep="\t")
    tab_new <- update_table(tab, tab_new, tags)
    if (verbose){
      cat(paste("updated table at", file), "\n")
    }
  }
  write.table(tab_new, file=file, sep="\t", row.names=F, quote=F)
}
