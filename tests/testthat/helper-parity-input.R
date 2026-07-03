# Deterministic inputs shared by parity tests.
parity_inputs <- function() {
  set.seed(42)
  env    <- 1:5
  trait  <- c(2.0, 3.5, 5.0, 6.5, 8.0)                              # clean linear
  trait2 <- c(1.0, 4.0, 2.0, 9.0, 3.0)                              # noisy
  gxe    <- matrix(c(2, 4, 6, 3, 6, 9, 1, 2, 3), nrow = 3, byrow = TRUE)  # 3 geno x 3 env

  # Replicated multi-environment design (3 environments x 3 replicates).
  env_r   <- rep(1:3, each = 3)
  trait_r <- c(2, 3, 4, 6, 7, 8, 3, 4, 5)
  group_r <- rep(c("g1", "g2", "g3"), times = 3)
  covar_r <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)

  # Two-timepoint series (n = 5) for developmental indices.
  time1 <- c(2, 3, 4, 5, 6)
  time2 <- c(3, 5, 4, 8, 7)

  # List-of-vectors (one per group) for the list-input index branches.
  trait_list <- list(g1 = c(2, 3, 4), g2 = c(6, 7, 8), g3 = c(3, 4, 5))

  # Misc scalars / vectors.
  baseline  <- c(2, 2, 2, 2, 2)
  binary5   <- c(0, 1, 0, 1, 0)
  trait_mat <- matrix(c(2, 4, 3, 6, 1, 2, 5, 7, 4, 3), ncol = 2, byrow = TRUE)  # 5 samples x 2 traits

  # Data frames for the data.frame-API indices.
  df_spi <- data.frame(env = env_r, y = trait_r)
  df_pr  <- data.frame(env = env_r, y = trait_r, pop = group_r)

  list(
    env = env, trait = trait, trait2 = trait2, gxe = gxe,
    env_r = env_r, trait_r = trait_r, group_r = group_r, covar_r = covar_r,
    trait_list = trait_list,
    time1 = time1, time2 = time2, baseline = baseline, binary5 = binary5,
    trait_mat = trait_mat, df_spi = df_spi, df_pr = df_pr
  )
}

# A named list of zero-arg thunks, one per index, each invoking `fns$<name>`
# (an environment or list of the calculate_* functions) on shared inputs.
# Thunks let the generator and the test evaluate calls individually and isolate
# any that error, so both exercise byte-identical calls.
#
# calculate_reaction_norm_non_linear is included: its return value depends only
# on the polynomial-fit coefficients, not on the dead `predict()` line that
# referenced an undefined `env` in the reference. The fixture generator supplies
# the global `env` the original monolithic pipeline relied on so the untouched
# reference runs and yields its genuine (equidistant) value.
parity_calls <- function(fns) {
  i <- parity_inputs()
  list(
    calculate_CVt                    = function() fns$calculate_CVt(i$trait2),
    calculate_reaction_norm_slope    = function() fns$calculate_reaction_norm_slope(i$trait, i$env),
    calculate_reaction_norm_non_linear = function() fns$calculate_reaction_norm_non_linear(i$trait2, degree = 2),
    calculate_D_slope                = function() fns$calculate_D_slope(i$trait2),
    calculate_RC                     = function() fns$calculate_RC(i$trait2),
    calculate_CVm                    = function() fns$calculate_CVm(i$trait_list),
    calculate_CVmd                   = function() fns$calculate_CVmd(i$trait_r, i$group_r),
    calculate_grand_plasticity       = function() fns$calculate_grand_plasticity(i$trait_r, i$env_r),
    calculate_PPF                    = function() fns$calculate_PPF(i$trait_r, i$env_r, i$covar_r),
    calculate_Phenotypic_Plasticity_Index = function() fns$calculate_Phenotypic_Plasticity_Index(i$trait2),
    calculate_PImd                   = function() fns$calculate_PImd(i$trait_r, i$env_r),
    calculate_PILSM                  = function() fns$calculate_PILSM(i$trait_r, i$env_r),
    calculate_RTR                    = function() fns$calculate_RTR(i$trait2, i$env),
    calculate_PIR                    = function() fns$calculate_PIR(i$trait_r, i$env_r),
    calculate_PSI                    = function() fns$calculate_PSI(i$trait2, i$env),
    calculate_RPI                    = function() fns$calculate_RPI(i$trait2, i$env),
    calculate_PQ                     = function() fns$calculate_PQ(i$trait2, i$env),
    calculate_PR                     = function() fns$calculate_PR(i$trait2, i$env),
    calculate_NRW                    = function() fns$calculate_NRW(i$trait_r, i$env_r, i$group_r),
    calculate_ESP                    = function() fns$calculate_ESP(i$trait2, i$env),
    calculate_general_PD             = function() fns$calculate_general_PD(i$trait2, i$env),
    calculate_FPI                    = function() fns$calculate_FPI(i$trait2, i$env),
    calculate_TPS                    = function() fns$calculate_TPS(i$trait_r, i$env_r, native_env = 1, transplanted_env = 3),
    calculate_DPI                    = function() fns$calculate_DPI(i$time1, i$time2),
    calculate_CEV                    = function() fns$calculate_CEV(i$trait2),
    calculate_PRI                    = function() fns$calculate_PRI(i$trait2, i$binary5),
    calculate_PFI                    = function() fns$calculate_PFI(i$trait2, i$baseline),
    calculate_SPI                    = function() fns$calculate_SPI(i$df_spi, "env", "y", env1 = 1, env2 = 3, reference_env = 2),
    calculate_APC                    = function() fns$calculate_APC(i$trait2, i$env),
    calculate_SI                     = function() fns$calculate_SI(i$trait2),
    calculate_RSI                    = function() fns$calculate_RSI(i$trait2),
    calculate_EVS                    = function() fns$calculate_EVS(i$trait2),
    calculate_MVPi                   = function() fns$calculate_MVPi(i$trait_mat),
    calculate_SPM                    = function() fns$calculate_SPM(i$df_spi, "env", "y", resident_env = 1, nonresident_env = 3),
    calculate_Plasticity_Ratio       = function() fns$calculate_Plasticity_Ratio(i$df_pr, "env", "y", "pop"),
    calculate_rdpi                   = function() fns$calculate_rdpi(i$trait2, i$env),
    calculate_ESPI                   = function() fns$calculate_ESPI(i$trait2, i$env),
    calculate_espiid                 = function() fns$calculate_espiid(i$trait2),
    calculate_plasticity             = function() fns$calculate_plasticity(i$trait2, i$env),
    calculate_finlay_wilkinson       = function() fns$calculate_finlay_wilkinson(i$gxe)
  )
}
