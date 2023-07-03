#' Transform Dataset to Pseudo Copula Variables
#'
#' This function calculates and applies a normal transformation on a per column basis, allowing for missing values.
#' The function calculates the forward and backwards transformation functions using an empirically smooth CDF to convert continuous data to pseudo gaussian variables.
#' It returns the inverse functions in order to back transform to the original data scale and the transformed data.
#'
#' @param dataset The dataset to transform select columns to gaussian variables
#' @param columns The specified columns to transform as column names
#' @param rownamevar The observation ID variable name to ensure the data are aligned when missing values are present.
#' @param iteration Either TRUE or FALSE to print the current column iteration to track progress or not.
#'
#' @return a list with transformed_df including the row IDs and all columns transformed to gaussian variables and qfunctions with all the inverse functions needed to backtransform to the original data scale
#' @export
#'
#' @examples
#' # determine which missing CpGs are non-normal at 0.01 threshold
#' Column_Normality_results <- TestNormalityofMissingCols(toy_DNAm, 0.01)
#' # extract column names into a vector
#' cols_below_threshold <- Column_Normality_results$Columns
#' # Transform the selected columns to Normal variables #
#' transformed_toy_DNAm <- TransformDataset(dataset = toy_DNAm,
#'                                           columns = cols_below_threshold,
#'                                           rownamevar = "PersonID", iteration = TRUE)
TransformDataset <- function(dataset, columns, rownamevar, iteration=TRUE){
  # Output list initialization
  col_output <- list()

  # Preparing the dataset for transformation
  row_name <- rownames(dataset)
  dataset_prep <- data.frame(dataset[, columns], ID = row_name)

  # Last column is the ID now
  last_val <- dim(dataset_prep)[2]
  # Define the last CpG column
  last_cpg <- last_val - 1

  # Loop through each column to apply transformations
  for(i in 1:last_cpg){
    col_name <- colnames(dataset_prep)[i]

    # Prepare the data: only non-missing values are processed.
    # The data frame now has ID in the first column and CpG in the second column
    dat_prep <- dataset_prep[!is.na(dataset_prep[, i]), c(last_val, i)]

    # Density estimation for non-missing CpG values, using kernel density estimation for smoothing
    # 2nd column used because the second column is the current CpG #
    density_x <- density(dat_prep[, 2], adjust = 0.5, from=min(dat_prep[, 2]),
                to=max(dat_prep[, 2]), n = 1000)

    # Transform the frequency into density function
    dapproxfun <- splinefun(x = density_x$x, y = density_x$y)
    dfun <- function(x) dapproxfun(x)
    support <- range(dat_prep[, 2])

    # Integrate the density function to generate the empirical CDF
    pfun_integrate_dfun_1 <- function(v) integrate(dfun, support[1], v,
                                                   subdivisions=2000, rel.tol =1e-15, stop.on.error = FALSE)$value
    pfun_integrate_dfun <- function(x) Vectorize(pfun_integrate_dfun_1)(x)

    # pfun is the empirical (smoothed) CDF
    pfun <- splinefun(x = dat_prep[, 2],
                  y = pfun_integrate_dfun(dat_prep[, 2]))

    # Transform to uniform space using the empirical smoothed CDF
    uniform_var <- pfun(dat_prep[, 2])

    # If some values are over 1, rescale and make sure to not have exact 1's or 0's
    if(sum(uniform_var > 1) > 0){
      uniform_var <- uniform_var / max(uniform_var)
      uniform_var <- ifelse(uniform_var == 1, 0.9999999, uniform_var)
    }
    uniform_var <- ifelse(uniform_var == 0, 0.0000001, uniform_var)

    # Compute the inverse of the empirical CDF (qfun)
    qfun <- splinefun(x = uniform_var, y = dat_prep[, 2])

    # Transform to standard normal space
    norm_var <- qnorm(uniform_var)
    # Handling infinite values by setting them to +/- 5
    norm_var <- ifelse(is.infinite(norm_var), sign(norm_var) * 5, norm_var)

    # Create dataframe with normal transformed variable and ID
    df_sm <- data.frame(norm_var, ID = dat_prep[, 1])
    colnames(df_sm) <- c(paste0(col_name, "_normal_var"), "ID")

    # Merge the original scale variable and normal_var
    df_sm2 <- merge(dataset_prep[, c("ID", col_name)], df_sm, by = "ID", all.x=TRUE)

    # Rename ID variable to original name
    colnames(df_sm2) <- c(paste0({{rownamevar}}), colnames(df_sm2)[2:3])

    # Prepare output for the current column
    col_output1 <- list(qfun, df_sm2)
    new_names <- c(paste0(col_name, "_qfun"), paste0(col_name, "_normal_var"))
    names(col_output1) <- new_names

    # Add the output to the main output list
    col_output[[i]] <- col_output1

    if(iteration == TRUE){
      print(i) # Print current iteration number for tracking progress
    }
  }

  # convert output to list of qfunctions and transformed dataframe
  qfunctions_output <- qfun_extract(col_output)

  # extract the transformed values
  transformed_df <- extract_transformed_values(col_output)

  output_return <- list(transformed_df, qfunctions_output)
  names(output_return) <- c("transformed_df", "qfunctions")
  return(output_return)
}
