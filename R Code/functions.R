# 1. imputedata----

#### Input:
#### - data = data to impute as matrix or dataframe
#### - m = Number of multiple imputations

#### Output:
#### List of m imputed dataframes. In each of these dataframes the missing values differ.

imputedata <- function(data, m = 20) {
  require("mice")
  
  # Imputate the data:
  impobj <- mice(data, m = m, printFlag = FALSE)
  dataimputed <- list()
  for(i in 1:m) {
    dataimputed[[i]] <- complete(impobj, action=i)
  }
  # Return the imputed data:
  return(dataimputed)
}

# 2. makeCVdiv----

#### Input:
#### - data = data for crossvalidation
#### - yname = Factor variable by which stratified cross-validation is done
#### - ncv = Number of folds drawn

#### Output:
#### Vector with the individual assignments for the dataframe

makeCVdiv <- function(data, yname="y", ncv=5) {
  cvid <- rep(NA, length = nrow(data))
  for(i in levels(data[, yname])) cvid[data[, yname] == i] <- sample(rep(1:ncv, length = sum(data[, yname] == i)))
  cvid
}