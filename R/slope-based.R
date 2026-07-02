#' Reaction norm slope
#'
#' Fits a simple linear regression of trait values on environmental values and
#' returns the slope, i.e. the linear rate of phenotypic change across
#' environments.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param environments Optional numeric vector of environment values, the same
#'   length as `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @return A single numeric value: the slope coefficient from
#'   `lm(trait_values ~ environments)`. Returns `NA` when fewer than two trait
#'   values are supplied.
#' @examples
#' calculate_reaction_norm_slope(c(2, 4, 6, 8))
#' @export
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

#' Non-linear reaction norm score
#'
#' Fits a raw polynomial regression of the given degree to trait values across
#' equidistant environments and summarizes non-linearity as the sum of the
#' absolute values of the non-intercept coefficients.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param degree Integer polynomial degree to fit. Defaults to `2`.
#' @return A single numeric value: the nonlinearity score. Returns `NA` with a
#'   warning when there are not enough data points for the requested degree.
#' @note The example below is illustrative only and is wrapped in `\dontrun{}`
#'   because the function body references an undefined `env` object when
#'   building `newdata` for `predict()`, rather than its own `environments`
#'   variable.
#' @examples
#' \dontrun{
#' calculate_reaction_norm_non_linear(c(2, 4, 9, 16, 25), degree = 2)
#' }
#' @export
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

#' Divergence slope (D_slope)
#'
#' Sorts trait values and computes the difference between the mean of the top
#' fraction and the mean of the bottom fraction, a robust measure of the spread
#' between the extremes of the trait distribution.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param lower_fraction Fraction of sorted values (from the bottom) defining the
#'   "low" segment. Defaults to `0.2`.
#' @param upper_fraction Fraction of sorted values (from the bottom) marking the
#'   start of the "high" segment. Defaults to `0.8`.
#' @return A single numeric value: `mean(upper_values) - mean(lower_values)`.
#' @examples
#' calculate_D_slope(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
#' @export
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

#' Response coefficient (RC)
#'
#' Sorts trait values, splits them into low and high segments, and computes the
#' ratio of the high-segment mean to the low-segment mean.
#'
#' @param trait_values Numeric vector of trait measurements.
#' @param lower_fraction Fraction of sorted values (from the bottom) defining the
#'   "low" segment. Defaults to `0.5`.
#' @param upper_fraction Fraction of sorted values (from the top) defining the
#'   "high" segment. Defaults to `0.5`.
#' @return A single numeric value: `mean(upper_values) / mean(lower_values)`.
#' @examples
#' calculate_RC(c(2, 4, 6, 8))
#' @export
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

#' Relative trait range (RTR)
#'
#' Compares mean trait values at the low and high ends of an environmental
#' gradient, normalized by the maximum absolute trait value.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Numeric vector of environment values, the same length as
#'   `trait_values`.
#' @param env_low Lower threshold for the environmental gradient. If `< 1`, it is
#'   treated as a quantile probability; otherwise as an absolute value. Defaults
#'   to `0.2`.
#' @param env_high Upper threshold for the environmental gradient (symmetric with
#'   `env_low`, i.e. the high cutoff is at the `1 - env_high` quantile when
#'   `env_high < 1`). Defaults to `0.2`.
#' @return A single numeric value: `(mean_high - mean_low) / max(abs(trait_values))`.
#' @examples
#' calculate_RTR(c(2, 4, 6, 8, 10), env_values = c(1, 2, 3, 4, 5))
#' @export
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

#' Plasticity index at relative growth rate maximum (PIR)
#'
#' Computes environment-level means of a trait, estimates (or accepts) relative
#' growth rates between successive environments, and expresses the trait range
#' relative to the mean at the environment of maximum relative growth rate.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Optional grouping vector of environment labels, the same
#'   length as `trait_values`. Defaults to equidistant indices `1, 2, ..., n`
#'   (each observation its own environment).
#' @param rgr_values Optional numeric vector of pre-computed relative growth
#'   rates, one per observation, averaged per environment. When `NULL`, relative
#'   growth rates are derived from successive environment means.
#' @return A single numeric value: `(max_mean - min_mean) / mean_at_max_rgr`.
#' @examples
#' calculate_PIR(c(2, 4, 6, 8))
#' @export
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

#' Best-fit polynomial plasticity score
#'
#' Fits polynomial regressions of `trait_values` on `env_values` for degrees `1`
#' through `max_degree`, selects the best-fitting model by AIC or BIC, and
#' summarizes plasticity as the sum of the absolute values of its non-intercept
#' coefficients.
#'
#' @param trait_values Numeric vector of trait measurements across environments.
#' @param env_values Optional numeric vector of environment values, the same
#'   length as `trait_values`. Defaults to equidistant indices `1, 2, ..., n`.
#' @param max_degree Maximum polynomial degree to consider. Defaults to `3`.
#' @param criterion Model selection criterion, `"BIC"` (default) or `"AIC"`.
#' @return A list with elements `best_degree`, `plasticity_score`, and
#'   `coefficients` (the non-intercept coefficients of the selected model).
#' @examples
#' calculate_plasticity(c(2.0, 3.5, 5.0, 6.5, 8.0), env_values = 1:5)
#' @export
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

#' Finlay-Wilkinson stability regression
#'
#' Regresses each genotype's trait values on an environmental index (by default
#' the centered environment means) to estimate Finlay-Wilkinson stability
#' (slope), intercept, fit statistics, and genotype effects. Optionally plots
#' the resulting regression lines.
#'
#' @param Y A genotype-by-environment numeric matrix of trait values (rows =
#'   genotypes, columns = environments).
#' @param genotype_ids Optional character or factor vector of genotype names, the
#'   same length as `nrow(Y)`. Defaults to the row names of `Y`, or
#'   `"G1", "G2", ...` if `Y` has none.
#' @param env_values Optional numeric vector of environmental covariate values,
#'   the same length as `ncol(Y)`. Defaults to the centered column means of `Y`.
#' @param plot Logical; if `TRUE`, plot the per-genotype regression lines.
#'   Defaults to `FALSE`.
#' @return A data frame with one row per genotype and columns `genotype`, `beta`
#'   (stability slope), `intercept`, `r2`, `rmse`, `n_env`, `G` (genotype effect),
#'   and `M` (grand mean).
#' @examples
#' Y <- matrix(c(2, 4, 7, 3, 6, 8, 1, 3, 4, 5, 9, 10), nrow = 4, byrow = TRUE)
#' calculate_finlay_wilkinson(Y)
#' @export
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
