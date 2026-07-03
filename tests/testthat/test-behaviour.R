test_that("constant trait yields zero variation for CVt", {
  x <- rep(5, 6)
  expect_equal(calculate_CVt(x), 0)
})

test_that("calculate_CVt rejects NA rather than silently tolerating it", {
  # calculate_CVt has no na.rm handling and explicitly stops on missing values;
  # this documents the actual (non-tolerant) current behaviour.
  x <- c(1, NA, 3, 5)
  expect_error(
    calculate_CVt(x),
    "missing values"
  )
})

test_that("calculate_ESP tolerates NA in the trait vector and computes the expected value", {
  # calculate_ESP uses na.rm = TRUE internally, so NA in one environment's
  # group should be dropped rather than propagating to NA/error.
  x <- c(1, NA, 3, 5)
  env <- c(1, 1, 2, 2)

  # Expected by hand:
  #  overall mean (na.rm) = mean(1, 3, 5) = 3
  #  env 1 mean (na.rm)   = mean(1)       = 1  -> deviation = (1 - 3) / 3 = -2/3
  #  env 2 mean           = mean(3, 5)    = 4  -> deviation = (4 - 3) / 3 =  1/3
  #  ESP = sum(abs(deviations)) = 2/3 + 1/3 = 1
  expect_silent(result <- calculate_ESP(x, env_values = env))
  expect_equal(result, 1)
})

test_that("rdpi_mean_calculation returns a long-format data frame of pairwise RDPI values", {
  df <- data.frame(env = rep(1:3, each = 2), t1 = c(1, 2, 2, 4, 3, 6))
  out <- rdpi_mean_calculation(df, trait_cols = "t1", factors = "env")

  # Observed shape: one row per unique pair of the 3 environment levels
  # (3 choose 2 = 3 rows), columns sp / env1 / env2 / rdpi.
  expect_s3_class(out, "data.frame")
  expect_equal(dim(out), c(3, 4))
  expect_equal(names(out), c("sp", "env1", "env2", "rdpi"))
  expect_equal(out$sp, rep("Single_Group", 3))
  expect_equal(out$env1, c("1", "1", "2"))
  expect_equal(out$env2, c("2", "3", "3"))

  # env means: env1 = mean(1,2) = 1.5, env2 = mean(2,4) = 3, env3 = mean(3,6) = 4.5
  # rdpi(env1,env2) = |1.5-3|/min(1.5,3)   = 1.0
  # rdpi(env1,env3) = |1.5-4.5|/min(1.5,4.5) = 2.0
  # rdpi(env2,env3) = |3-4.5|/min(3,4.5)   = 0.5
  expect_equal(out$rdpi, c(1.0, 2.0, 0.5))
})

test_that("calculate_reaction_norm_non_linear returns the nonlinearity score", {
  # Regression for the fixed `env`-not-found bug: the function used to reference
  # an undefined `env` object in dead code and errored for every input. The
  # score is the sum of |coefficients| of the raw polynomial fit, excluding the
  # intercept. Oracle computed independently from a plain lm().
  y <- c(1, 4, 2, 9, 3)
  coefs <- coef(lm(y ~ poly(seq_along(y), 2, raw = TRUE)))
  expect_equal(calculate_reaction_norm_non_linear(y, degree = 2), sum(abs(coefs[-1])))
})

test_that("calculate_reaction_norm_non_linear honours non-equidistant environments", {
  # Supplying environment spacing changes the polynomial fit (and hence the
  # score); NULL falls back to equidistant seq_along().
  y   <- c(1, 4, 2, 9, 3)
  env <- c(1, 2, 4, 8, 16)
  coefs <- coef(lm(y ~ poly(env, 2, raw = TRUE)))
  expect_equal(calculate_reaction_norm_non_linear(y, environments = env, degree = 2),
               sum(abs(coefs[-1])))
  # Explicit equidistant spacing matches the default.
  expect_equal(calculate_reaction_norm_non_linear(y, environments = seq_along(y)),
               calculate_reaction_norm_non_linear(y))
  # Length mismatch is rejected.
  expect_error(calculate_reaction_norm_non_linear(y, environments = c(1, 2, 3)),
               "same length")
})

test_that("calculate_CVm accepts a vector plus env_values", {
  # Regression for the fixed `group_labels`-not-found bug: the non-list branch
  # referenced an undefined `group_labels` object instead of `env_values`.
  trait <- c(2, 3, 4, 6, 7, 8, 3, 4, 5)
  env   <- rep(1:3, each = 3)
  m <- tapply(trait, env, mean)          # group means: 3, 7, 4
  expect_equal(calculate_CVm(trait, env_values = env), sd(m) / mean(m))
})

test_that("calculate_PPF works without covariates", {
  # Regression for the fixed emmeans "reference grid" bug: the no-covariate
  # branch fitted a model whose term name did not match the emmeans spec.
  # With a balanced one-way design the LSMeans equal the group means, so the
  # expected value can be derived without emmeans.
  trait <- c(2, 3, 4, 6, 7, 8, 3, 4, 5)
  env   <- rep(1:3, each = 3)
  m <- tapply(trait, env, mean)          # 3, 7, 4
  pairs <- combn(names(m), 2, simplify = FALSE)
  # PPF per pair = 100 * |LSM1 - LSM2| / LSM1, averaged over all pairs.
  expected <- mean(vapply(pairs, function(p) 100 * abs(m[[p[1]]] - m[[p[2]]]) / m[[p[1]]], numeric(1)))
  expect_equal(calculate_PPF(trait, env), expected)
})
