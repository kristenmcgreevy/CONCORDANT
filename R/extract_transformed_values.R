#' Extract the Normal Transformed Values from the Transformation Procedure
#'
#' @param input_list
#'
#' @return dataframe with row ID column and Normal transformed variables as other columns
#'
#'
#' @examples
extract_transformed_values <- function(input_list){

  transformed_df <- data.frame(ID = input_list[[1]][[2]][, 1])
  colnames(transformed_df) <- colnames(input_list[[1]][[2]])[1]

  for (i in seq_len(length(input_list))) {
    transformed_var <- input_list[[i]][[2]][, 3]
    colname <- colnames(input_list[[i]][[2]])[3]
    transformed_df[colname] <- transformed_var
  }

  return(transformed_df)
}
