# CONCORDANT
COpula Normalization and COntinuous Data ANalysis Toolkit

These R functions and source code calculate functions to transform variables to pseudo gaussian variables using copula transformations. It allows researchers to change the distribution of their data to better align with common statistical tools necessitating normally distributed variables. This package accompanies a manuscript that demonstrates its use with DNAm data, which is commonly skewed or bimodal and restricted to values between 0 and 1. 

The `ExampleRScript.R` holds code to run through how to use these functions on a small, simulated DNAm dataframe with missing values. 

#### `TestNormalityofMissingCols` ####
This function finds columns with at least 1 missing value and evaluates whether it is normally distributed or not using Shapiro Wilks test at your specified p-value threshold. 

To use this function, supply 2 arguments: 

`dataset`: dataset with variables in columns and observations in rows

`p_value_threshold`: either `'any', 0.05, 0.01 or 0.001`. When set to "any", the function will return the names of all columns with missing values, regardless of the p-value from the Shapiro-Wilks test.

The function returns a dataframe two columns: 

`Columns` contains the names of the columns with p-values below the specified threshold, and 

`P_values` contains the corresponding p-values from the Shapiro-Wilk tests.



#### `TransformDataset` ####
Transform your dataset to pseudo gaussian variables. This function calculates and applies a normal transformation on a per column basis, allowing for missing values. The function calculates the forward and backwards transformation functions using an empirically smooth CDF to convert continuous data to pseudo gaussian variables.

To use this function, supply 4 arguments:

`dataset`: The dataset to transform select columns to gaussian variables

`columns`: The specified columns to transform (as column names)

`rownamevar`: The observation ID variable to ensure the data are aligned when missing values are present.

`iteration`: Either `TRUE` or `FALSE` to print the current iteration to track progress.

The function returns a 2 item list:

`transformed_df`: has the row IDs and columns transformed to gaussian variables. Note that columns transformed have `_normal_var` added to column names

`qfunctions`: the inverse functions needed to backtransform the data


#### `BackTransformDataset` #### 
Back transform your normal transformed data to the original data scale. This function applies the inverse CDF (q-function) of each column to back-transform. It only applies the backtransformation to columns with supplied q-functions.
To use this function, supply 2 arguments:
`dataset`: The dataset on the gaussian scale to transform select columns back to variables on the original data scale
`qfunctions`: A list of the inverse functions to apply column-wise
The function returns a dataframe with original variable names and variables on the original datascale.

