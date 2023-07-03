#' Test Normality of Missing Columns
#'
#' Evaluate all columns in the dataset with at least 1 missing value for being normally distributed or not. Shapiro Wilks p-value is used where p-values under the threshold are deemed non-normally distributed and p-values above the threshold as deemed relatively normally distributed.
#'
#' This function is intended to be used for the later functions to know which columns to transform to pseudo gaussian variables.
#'
#' @param dataset dataset with variables in columns and observations in rows
#' @param p_value_threshold either 'any', 0.05, 0.01 or 0.001. This threshold corresponds to the Shapiro Wilks p-value threshold for determining if the column is normally distributed or not. The function also handles the case where the threshold is set to "any", meaning that it will return the names of all columns with missing values, regardless of the p-value from the Shapiro-Wilks test.
#'
#' @return dataframe with two columns: "Columns" contains the names of the columns with p-values below the specified threshold, and "P_values" contains the corresponding p-values from the Shapiro-Wilk tests.
#' @export
#'
#' @examples
#' # determine which missing CpGs are non-normal at 0.01 threshold
#' Column_Normality_results <- TestNormalityofMissingCols(toy_DNAm, 0.01)
TestNormalityofMissingCols <- function(dataset, p_value_threshold = 0.001) {

  if(!(p_value_threshold %in% c(0.05, 0.01, 0.001, "any"))) {
    stop("p_value_threshold must be either 'any', 0.05, 0.01 or 0.001")
  }

  # Identify columns with missing values
  cols_with_na <- colnames(dataset)[colSums(is.na(dataset)) > 0]

  # List to store columns with p-values below the specified threshold
  cols_below_threshold <- vector()
  p_values_below_threshold <- vector()

  # Loop through each column with missing values
  for(col_name in cols_with_na) {

    # Exclude missing values
    col_values <- dataset[[col_name]][!is.na(dataset[[col_name]])]

    # Only calculate Shapiro-Wilk test for column with at least 3 unique values
    if(length(unique(col_values)) >= 3) {
      # Calculate the Shapiro-Wilk test
      shapiro_test <- shapiro.test(col_values)

      # If p-value is below the specified threshold, add column name and p-value to lists
      if((shapiro_test$p.value < p_value_threshold) | (p_value_threshold == "any")) {
        cols_below_threshold <- c(cols_below_threshold, col_name)
        p_values_below_threshold <- c(p_values_below_threshold, shapiro_test$p.value)
      }
    }
  }

  output <- data.frame(Columns = cols_below_threshold, P_values = p_values_below_threshold)
  return(output)
}
