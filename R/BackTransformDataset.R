#' Back Transform Normal Transformed Data to the Original Scale
#'
#' Function for back-transforming data to its original scale. This assumes the data is in standard normal scale. It applies the inverse CDF (q-function) of each column to back-transform. It only applies the backtransformation to columns with supplied q-functions.
#'
#' @param dataset The dataset to transform select gaussian columns to variables on the original scale
#' @param qfunctions A list of the inverse functions to apply column-wise
#'
#' @return
#' @export
#'
#' @examples
#' # determine which missing CpGs are non-normal at 0.01 threshold
#' Column_Normality_results <- TestNormalityofMissingCols(toy_DNAm, 0.01)
# extract column names into a vector
#' cols_below_threshold <- Column_Normality_results$Columns
#'
# Transform the selected columns to Normal variables #
#' transformed_toy_DNAm <- TransformDataset(dataset = toy_DNAm,
#'                                           columns = cols_below_threshold,
#'                                           rownamevar = "PersonID", iteration = TRUE)
#'
#' # extract the qfunctions
#' missingqfunctions <- transformed_toy_DNAm$qfunctions
#'
#' # extract the transformed dataframe
#' missingnormal_toyDNAm <- transformed_toy_DNAm$transformed_df
#'
#' # merge with untransformed data
#' match1 <- match(as.numeric(missingnormal_toyDNAm$PersonID), toy_DNAm$PersonID)
#'
#' toy_DNAm_4imputation <- data.frame(missingnormal_toyDNAm,
#'                                    toy_DNAm[match1, colnames(toy_DNAm) %!in%
#'                                               c("PersonID", cols_below_threshold)])
#'
#' library(impute) # needs to be a matrix, and needs to have ID variable removed or set to numeric value
#' toy_DNAm_4imputation$PersonID <- as.numeric(toy_DNAm_4imputation$PersonID)
#' toy_DNAm_4imputation_unlist <- as.numeric(unlist(toy_DNAm_4imputation))
#' MissNormal_missing_df2 <- matrix(toy_DNAm_4imputation_unlist,
#'                                  nrow = nrow(toy_DNAm_4imputation),
#'                                  ncol = ncol(toy_DNAm_4imputation))
#' colnames(MissNormal_missing_df2) <- colnames(toy_DNAm_4imputation)
#' rownames(MissNormal_missing_df2) <- toy_DNAm_4imputation$PersonID
#'
#' knn_DNAm_impute <- impute.knn(MissNormal_missing_df2, k = 50)
#'
#' # get imputed dataframe (on Normal scale) #
#' imputed_toyDNAm <- knn_DNAm_impute$data
#' # Finally, back transform the imputed dataset to the original CpG beta values
#' backtrans_toyDNAm <- BackTransformDataset(imputed_dataset = imputed_toyDNAm,
#'                                                qfunctions = missingqfunctions)
BackTransformDataset <- function(dataset, qfunctions){

  # Re-order the columns of the imputed data to match the order of qfunctions
  newqfuncname <- gsub("*_qfun", "_normal_var", names(qfunctions))
  col_neworder <- match(newqfuncname, colnames(dataset))
  impute_df_matched <- dataset[, col_neworder]

  # Identify the remaining columns that were imputed but don't need to be transformed
  # These will be combined with the back-transformed columns later
  rest_imputed <- colnames(dataset)[colnames(dataset) %!in% newqfuncname]
  impute_df_rest <- dataset[, rest_imputed]


  # Start back transformation process
  col_num <- length(qfunctions)
  row_num <- nrow(impute_df_matched)

  # Apply the standard normal CDF (pnorm) to the imputed data
  # This transforms it to a uniform distribution
  missnormal_rep_unif <- pnorm(impute_df_matched)

  # Pre-allocate a dataframe for the back-transformed data
  backtransform_missnormal_rep <- data.frame(matrix(NA, ncol = col_num, nrow = row_num))

  # For each column, apply the corresponding inverse CDF (q-function)
  for(i in 1:col_num){
    backtransform_missnormal_rep[, i] <- qfunctions[[i]](missnormal_rep_unif[, i])
  }

  # Restore the original column names
  new_colnames <- gsub("*_normal_var", "\\1", colnames(missnormal_rep_unif))
  colnames(backtransform_missnormal_rep) <- new_colnames

  # Combine the back-transformed columns with the remaining imputed columns
  impute_df_merged <- cbind(backtransform_missnormal_rep, impute_df_rest)

  return(impute_df_merged)
}
