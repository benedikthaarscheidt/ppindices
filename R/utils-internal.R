impute = function(mt, mode = "median") {
  if (!is.matrix(mt) && !is.data.frame(mt)) {
    stop("Input must be a matrix or data frame")
  }

  ismt = is.matrix(mt)
  if (ismt) {
    mt = as.data.frame(mt)
  }

  if (mode == "median") {
    for (i in 1:ncol(mt)) {
      if (!is.factor(mt[, i])) {
        mt[is.na(mt[, i]), i] = median(mt[, i], na.rm = TRUE)
      }
    }
  } else if (mode == "mean") {
    for (i in 1:ncol(mt)) {
      if (!is.factor(mt[, i])) {
        mt[is.na(mt[, i]), i] = mean(mt[, i], na.rm = TRUE)
      }
    }
  } else if (mode == "rpart") {
    if (nrow(mt) > 1000) {
      warning("Imputation using rpart may take a while for large datasets.")
    }
    for (i in 1:ncol(mt)) {
      midx = which(is.na(mt[, i]))
      if (length(midx) == 0) next
      idx = which(!is.na(mt[, i]))
      colname = colnames(mt)[i]
      frm = as.formula(paste(colname, "~ ."))
      mod = rpart(frm, data = mt[idx, ], method = ifelse(is.factor(mt[, i]), "class", "anova"))
      vals = predict(mod, newdata = mt[midx, ], type = ifelse(is.factor(mt[, i]), "class", "vector"))
      mt[midx, i] = vals
    }
  } else if (mode == "knn") {
    if (!requireNamespace("VIM", quietly = TRUE)) {
      choice = utils::menu(c("Yes", "No"), title = "Package 'VIM' is required for k-NN imputation. Would you like to install it?")
      if (choice == 1) {
        install.packages("VIM")
      } else {
        stop("k-NN imputation requires the 'VIM' package. Please install it to proceed.")
      }
    }
    mt = VIM::kNN(mt, k = 5, imp_var = FALSE)
  } else {
    stop(paste("mode", mode, "not yet implemented"))
  }

  if (ismt) {
    mt = as.matrix(mt)
  }

  return(mt)
}

################################

combine_factors = function(dataframe, factors = NULL, factors_not_in_dataframe = NULL) {
  # Combine internal and external factors into a single dataframe
  if (!is.null(factors_not_in_dataframe)) {
    # Ensure the lengths match
    if (length(factors_not_in_dataframe[[1]]) != nrow(dataframe)) {
      stop("The length of external factors must match the number of rows in the dataframe.")
    }
    # Create a data frame for external factors
    external_factors_df = as.data.frame(factors_not_in_dataframe)

    # If there are internal factors, combine them with external factors
    if (!is.null(factors)) {
      factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
      combined_factors_df = cbind(dataframe[factors], external_factors_df)
    } else {
      combined_factors_df = external_factors_df
    }

    # Create a combined factor interaction
    dataframe$Combined_Factors = interaction(combined_factors_df, drop = TRUE)
  } else if (!is.null(factors)) {
    # If only internal factors are provided
    factors = if (is.numeric(factors)) names(dataframe)[factors] else factors
    dataframe$Combined_Factors = interaction(dataframe[factors], drop = TRUE)
  } else {
    stop("You must provide either internal factors, external factors, or both.")
  }

  # Ensure Combined_Factors is treated as a factor
  dataframe$Combined_Factors = as.factor(dataframe$Combined_Factors)
  levels=levels(dataframe$Combined_Factors)
  print(levels)
  return(dataframe)
}
