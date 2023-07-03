###################################################
### Source Functions for Copula Transform Paper ###
### Updated by Kristen July 3 2023              ###
###################################################

# Not in from dplyr 
`%!in%` = Negate(`%in%`)

# This function tests for normality (shapiro wilks) of the missing columns 
# based on the specified p-value threshold. It returns the column names that do
# not adhere to the threshold and their corresponding p-value 
# This version of the function now returns a dataframe with two columns: 
# "Columns" contains the names of the columns with p-values below the specified threshold, and 
# "P_values" contains the corresponding p-values from the Shapiro-Wilk tests.
# The function also handles the case where the threshold is set to "any",  
# meaning that it will return the names of all columns with missing values, 
# regardless of the p-value from the Shapiro-Wilk test.
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
  # names(output) <- c("Columns", "P_values")
  return(output)
}


# The function transform_dataset applies a normal transformation 
# to a dataset on a per column basis, handling missing values and creating an empirical 
# cumulative distribution function (CDF). 
# It outputs the inverse CDF for each column and the transformed data.
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


# Function for back-transforming imputed data to its original scale
# Assumes the imputed data is in standard normal scale (mean 0, sd 1)
# Applies the inverse CDF (q-function) of each column to back-transform
# Only applies to columns for which a q-function is provided (missing normal columns)
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



# convert my list outputs with normalized values into one large list  
# this is internal function to be called / used by other functions above 
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

# get all the inverse functions #
# this is internal function to be called / used by other functions above 
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
