calculate_reaction_norm_slope = function(trait_values, environments=NULL) {
  if (length(trait_values) < 2) {
    return(NA)
  }

  if(is.null(environments)){
    environments = seq_along(trait_values)
  }

  lm_fit = lm(trait_values ~ environments)

  return(unname(coef(lm_fit)["environments"]))
}

################################

calculate_reaction_norm_non_linear <- function(trait_values, degree = 2) {
  environments <- seq_along(trait_values)

  # Ensure there are enough data points for the model
  if (length(trait_values) < degree + 1) {
    warning("Not enough data points for the specified degree; returning NA.")
    return(NA)
  }

  # Fit the polynomial model using raw polynomial terms
  model <- tryCatch({
    lm(trait_values ~ poly(environments, degree, raw = TRUE))
  }, error = function(e) {
    warning("Model fitting failed; returning NA.")
    return(NULL)
  })

  if (is.null(model)) return(NA)
  fitted_vals <- predict(model, newdata = data.frame(env = env))
  range_val    <- max(fitted_vals) - min(fitted_vals)
  # Extract coefficients (first coefficient is the intercept)
  coefs <- coef(model)

  # Compute nonlinearity score as the sum of the absolute values of the coefficients for terms > 0
  nonlinearity_score <- sum(abs(coefs[-1]))

  return(nonlinearity_score)
}

################################

# NOTE: if the resource availability is an actual measurement of a metabolite being used by the
# plant, then grouping plants into high vs low resource availability should be done by clustering.
calculate_D_slope = function(trait_values, lower_fraction = 0.2, upper_fraction = 0.8) {
  # Ensure trait values are sorted and valid
  sorted_values = sort(trait_values, na.last = TRUE)

  # Calculate boundaries
  lower_boundary_index = floor(length(sorted_values) * lower_fraction)
  upper_boundary_index = ceiling(length(sorted_values) * upper_fraction)

  # Extract lower and upper segments
  lower_values = sorted_values[1:lower_boundary_index]
  upper_values = sorted_values[upper_boundary_index:length(sorted_values)]

  # Calculate D slope
  D_slope = mean(upper_values, na.rm = TRUE) - mean(lower_values, na.rm = TRUE)

  return(D_slope)
}

################################

calculate_RC = function(trait_values, lower_fraction = 0.5, upper_fraction = 0.5) {
  # Ensure trait values are sorted
  sorted_values = sort(trait_values, na.last = TRUE)

  # Calculate indices for splitting
  lower_boundary_index = floor(length(sorted_values) * lower_fraction)
  upper_boundary_index = ceiling(length(sorted_values) * (1 - upper_fraction))

  # Extract lower and upper segments
  lower_values = sorted_values[1:lower_boundary_index]
  upper_values = sorted_values[(upper_boundary_index + 1):length(sorted_values)]

  # Calculate means
  mean_low = mean(lower_values, na.rm = TRUE)
  mean_high = mean(upper_values, na.rm = TRUE)

  # Calculate Response Coefficient
  RC = mean_high / mean_low

  return(RC)
}

################################

calculate_RTR = function(trait_values, env_values, env_low = 0.2, env_high = 0.2) {
  # Validate input lengths
  if (length(trait_values) != length(env_values)) {
    stop("`trait_values` and `env_values` must have the same length.")
  }


  if (env_low < 1) {
    low_threshold = quantile(env_values, probs = env_low, na.rm = TRUE)
  } else {
    low_threshold = env_low
  }

  if (env_high < 1) {
    high_threshold = quantile(env_values, probs = 1 - env_high, na.rm = TRUE)
  } else {
    high_threshold = env_high
  }

  # Subset data based on thresholds
  low_idx = which(env_values <= low_threshold)
  high_idx = which(env_values >= high_threshold)

  # Calculate mean trait values for each end of the gradient
  mean_low = mean(trait_values[low_idx], na.rm = TRUE)
  mean_high = mean(trait_values[high_idx], na.rm = TRUE)

  # Calculate the RTR value
  RTR_value = (mean_high - mean_low) / max(abs(trait_values), na.rm = TRUE)

  return(RTR_value)
}

################################

calculate_PIR = function(trait_values, env_values = NULL, rgr_values = NULL) {

  if (is.null(env_values)) {
    env_values = factor(seq_along(trait_values))
  } else {
    env_values = as.factor(env_values)
  }


  means = tapply(trait_values, env_values, mean, na.rm = TRUE)


  if (is.null(means) || length(means) != length(levels(env_values))) {
    stop("Mismatch between environments and trait values. Ensure proper input alignment.")
  }


  max_mean = max(means, na.rm = TRUE)
  min_mean = min(means, na.rm = TRUE)


  if (is.null(rgr_values)) {
    rgr_values = numeric(length(means))
    rgr_values[2:length(means)] = diff(means) / means[-length(means)]
    rgr_values[1] = NA
  } else {

    rgr_values = tapply(rgr_values, env_values, mean, na.rm = TRUE)
  }


  if (length(rgr_values) != length(means)) {
    stop("Mismatch between calculated RGR values and environment means.")
  }


  max_rgr_index = which.max(rgr_values)


  mean_at_max_rgr = unname(means[max_rgr_index])


  PIR_value = (max_mean - min_mean) / mean_at_max_rgr

  return(PIR_value)
}

################################

calculate_plasticity = function(trait_values, env_values = NULL, max_degree = 3, criterion = "BIC") {
  # Input validation
  if (!is.numeric(trait_values)) {
    stop("trait_values must be a numeric vector")
  }

  n_values = length(trait_values)

  # If no environmental values are provided, create equidistant values
  if (is.null(env_values)) {
    env_values = seq_len(n_values)
  }

  # Ensure environmental values are numeric
  if (!is.numeric(env_values)) {
    stop("env_values must be numeric")
  }

  # Check for length mismatch
  if (length(trait_values) != length(env_values)) {
    stop("trait_values and env_values must have the same length")
  }

  # Store best model information
  best_degree = 1
  best_criterion_value = Inf
  best_model = NULL

  # Try polynomial degrees from 1 to max_degree
  for (degree in 1:max_degree) {
    formula = as.formula(paste("trait_values ~ poly(env_values, ", degree, ", raw = TRUE)", sep = ""))
    model = lm(formula)

    # Use AIC or BIC for model selection
    model_criterion = ifelse(criterion == "AIC", AIC(model), BIC(model))

    # Penalize higher-degree models if difference is small (<2)
    if ((model_criterion) < best_criterion_value) {
      best_criterion_value = model_criterion
      best_degree = degree
      best_model = model
    }
  }


  coefficients = coef(best_model)[-1]

  # Compute plasticity score as sum of absolute values of coefficients
  plasticity_score = sum(abs(coefficients))

  return(list(
    best_degree = best_degree,
    plasticity_score = plasticity_score,
    coefficients = coefficients
  ))
}

################################

calculate_finlay_wilkinson <- function(Y, genotype_ids=NULL, env_values=NULL, plot=FALSE) {
  Y <- as.matrix(Y)
  if (is.null(genotype_ids)) {
    gnames <- rownames(Y)
    if (is.null(gnames)) gnames <- paste0("G", seq_len(nrow(Y)))
  } else {
    if (length(genotype_ids) != nrow(Y)) stop("genotype_ids length must match nrow(Y)")
    gnames <- as.character(genotype_ids)
  }
  M <- mean(Y, na.rm=TRUE)
  if (is.null(env_values)) {
    X <- colMeans(Y, na.rm=TRUE) - M
    xlab <- "Environment index E_j (centered)"
  } else {
    if (length(env_values) != ncol(Y)) stop("env_values length must match ncol(Y)")
    X <- as.numeric(env_values) - mean(as.numeric(env_values), na.rm=TRUE)
    xlab <- "Covariate (centered)"
  }
  G_eff <- rowMeans(Y, na.rm=TRUE) - M
  res_list <- lapply(seq_len(nrow(Y)), function(i) {
    y <- as.numeric(Y[i, ])
    ok <- is.finite(y) & is.finite(X)
    x <- X[ok]; yy <- y[ok]
    if (length(yy) < 2 || var(x, na.rm=TRUE) == 0) {
      return(data.frame(genotype=gnames[i], beta=NA_real_, intercept=NA_real_, r2=NA_real_, rmse=NA_real_, n_env=length(yy), G=G_eff[i], M=M))
    }
    m <- lm(yy ~ x)
    pr <- predict(m)
    data.frame(genotype=gnames[i], beta=coef(m)[2], intercept=coef(m)[1], r2=summary(m)$r.squared, rmse=sqrt(mean((yy - pr)^2)), n_env=length(yy), G=G_eff[i], M=M)
  })
  res <- do.call(rbind, res_list); rownames(res) <- res$genotype
  if (plot) {
    plot(NA, xlim=range(X, na.rm=TRUE), ylim=range(Y, na.rm=TRUE), xlab=xlab, ylab="Trait Y_ij", main="Finlay–Wilkinson")
    cols <- setNames(seq_len(nrow(Y)), gnames)
    for (i in seq_len(nrow(Y))) {
      y <- as.numeric(Y[i, ])
      ok <- is.finite(y) & is.finite(X)
      points(X[ok], y[ok], pch=19, col=cols[gnames[i]])
      if (!is.na(res[i, "beta"])) abline(a=res[i, "intercept"], b=res[i, "beta"], col=cols[gnames[i]], lwd=2)
    }
    legend("topleft", legend=gnames, col=seq_len(nrow(Y)), pch=19, lwd=2, bty="n", title="Genotype")
  }
  res[, c("genotype","beta","intercept","r2","rmse","n_env","G","M")]
}
