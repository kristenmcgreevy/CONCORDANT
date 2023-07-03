#' Extract Inverse Functions from Transformation Output
#'
#' @param list_vals
#'
#' @return a list with each element the qfunction to backtransform the normal column to the original data scale.
#'
#' @examples
qfun_extract <- function(list_vals){

  qfunc_names <- vector()
  return_list <- list()
  res_1_list <- list()

  for(j in 1:length(list_vals)){
    # res_1 this has q function and transformed data for jth cpg
    res_1 <- list_vals[[j]]

    # this is just q function, which is first list element
    res_1_list[[j]] <- res_1[[1]]
    # get name of the qfunction
    qfunc_names[j] <- names(list_vals[[j]][1])

    # jth item is inverse function
    return_list[[j]] <- res_1_list[[j]]
  }
  # rename the list
  names(return_list) <- qfunc_names
  return(return_list)
}
