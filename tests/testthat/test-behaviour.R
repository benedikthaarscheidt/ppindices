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
