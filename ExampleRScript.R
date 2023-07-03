################################################################################
### Example R Script for Copula Transforms on DNAm data                      ###
### Run this code after downloading all files into the same directory        ###
################################################################################




################################################################################
###                          GETTING STARTED                                 ###  
################################################################################

# make sure these packages are installed and loaded #
library(dplyr)
library(impute)

# This source call will load in 6 functions into your global environment:
# TestNormalityofMissingCols, TransformDataset, BackTransform_ImputedData, 
# and internal functions: extract_transformed_values, qfun_extract, and %!in%
source("CopulaTransformSourceFunctions.R")

# if you are interested in working with a Toy DNAm dataset, load the data:
toy_DNAm <- readRDS("ToyDNAmData.rds")




################################################################################
###                        ABOUT THE FUNCTIONS                               ###  
################################################################################

###   TestNormalityofMissingCols   ###
# This function finds columns with at least 1 missing value and evaluates whether 
# it is normally distributed or not using Shapiro Wilks test at your specified
# p-value threshold. 
# To use this function, supply 2 arguments: 
# dataset: dataset with variables in columns and observations in rows
# p_value_threshold: either 'any', 0.05, 0.01 or 0.001. 
# When set to "any", the function will return the names of all columns with 
# missing values, regardless of the p-value from the Shapiro-Wilks test.
# The function returns a dataframe two columns: "Columns" contains the names of 
# the columns with p-values below the specified threshold, and "P_values" 
# contains the corresponding p-values from the Shapiro-Wilk tests.



###   TransformDataset   ###
# Transform your dataset to pseudo gaussian variables. This function calculates 
# and applies a normal transformation on a per column basis, allowing for 
# missing values. The function calculates the forward and backwards 
# transformation functions using an empirically smooth CDF to convert continuous 
# data to pseudo gaussian variables.
# To use this function, supply 4 arguments:
# dataset: The dataset to transform select columns to gaussian variables
# columns: The specified columns to transform (as column names)
# rownamevar: The observation ID variable to ensure the data are aligned when 
#             missing values are present.
# iteration: Either TRUE or FALSE to print the current iteration to track progress.
# The function returns a 2 item list:
# transformed_df: has the row IDs and columns transformed to gaussian variables. 
#                 Note that columns transformed have _normal_var added to column names
# qfunctions: the inverse functions needed to backtransform the data



###   BackTransformDataset   ###
# Back transform your normal transformed data to the original data scale. This
# function applies the inverse CDF (q-function) of each column to back-transform. 
# It only applies the backtransformation to columns with supplied q-functions.
# To use this function, supply 2 arguments:
# dataset: The dataset on the gaussian scale to transform select columns back 
#          to variables on the original data scale
# qfunctions: A list of the inverse functions to apply column-wise
# The function returns a dataframe with original variable names and variables on 
# the original datascale




################################################################################
###                            EXAMPLE CODE                                 ###
################################################################################

# determine which missing CpGs are non-normal at 0.01 threshold
Column_Normality_results <- TestNormalityofMissingCols(toy_DNAm, 0.01)
# extract column names into a vector
cols_below_threshold <- Column_Normality_results$Columns

# Transform the selected columns to Normal variables # 
transformed_toy_DNAm <- TransformDataset(dataset = toy_DNAm, 
                                          columns = cols_below_threshold, 
                                          rownamevar = "PersonID", iteration = TRUE)

# extract the qfunctions 
missingqfunctions <- transformed_toy_DNAm$qfunctions

# extract the transformed dataframe 
missingnormal_toyDNAm <- transformed_toy_DNAm$transformed_df

# merge with untransformed data
match1 <- match(as.numeric(missingnormal_toyDNAm$PersonID), toy_DNAm$PersonID)

toy_DNAm_4imputation <- data.frame(missingnormal_toyDNAm, 
                                   toy_DNAm[match1, colnames(toy_DNAm) %!in% 
                                              c("PersonID", cols_below_threshold)])


# Imputation # 
# needs to be a matrix for imputation, and needs to have ID variable removed or set to numeric value
toy_DNAm_4imputation$PersonID <- as.numeric(toy_DNAm_4imputation$PersonID)
toy_DNAm_4imputation_unlist <- as.numeric(unlist(toy_DNAm_4imputation))
MissNormal_missing_df <- matrix(toy_DNAm_4imputation_unlist, 
                                 nrow = nrow(toy_DNAm_4imputation), 
                                 ncol = ncol(toy_DNAm_4imputation))
colnames(MissNormal_missing_df) <- colnames(toy_DNAm_4imputation)
rownames(MissNormal_missing_df) <- toy_DNAm_4imputation$PersonID

knn_DNAm_impute <- impute.knn(MissNormal_missing_df, k = 50)

# get imputed dataframe (on Normal scale) #
imputed_toyDNAm <- knn_DNAm_impute$data

# backtransform data to original DNAm beta scale 
backtrans_toyDNAm <- BackTransformDataset(dataset = imputed_toyDNAm, 
                                               qfunctions = missingqfunctions)
# get columns in original order 
backtrans_toyDNAm2 <- backtrans_toyDNAm[, colnames(toy_DNAm)]

# save dataset for downstream analyses
saveRDS(backtrans_toyDNAm2, "BacktransformedDNAmData.rds")


