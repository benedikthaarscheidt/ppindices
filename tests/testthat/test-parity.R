test_that("packaged indices reproduce reference outputs bit-for-bit", {
  expected <- readRDS(test_path("fixtures", "parity.rds"))
  i <- parity_inputs()
  got <- suppressWarnings(list(
    calculate_CVt                 = calculate_CVt(i$trait2),
    calculate_reaction_norm_slope = calculate_reaction_norm_slope(i$trait, i$env),
    calculate_D_slope             = calculate_D_slope(i$trait2),
    calculate_RC                  = calculate_RC(i$trait2),
    calculate_CEV                 = calculate_CEV(i$trait2),
    calculate_PSI                 = calculate_PSI(i$trait2, i$env),
    calculate_PQ                  = calculate_PQ(i$trait2, i$env),
    calculate_PR                  = calculate_PR(i$trait2, i$env),
    calculate_ESP                 = calculate_ESP(i$trait2, i$env),
    calculate_SI                  = calculate_SI(i$trait2),
    calculate_RSI                 = calculate_RSI(i$trait2),
    calculate_EVS                 = calculate_EVS(i$trait2),
    calculate_rdpi                = calculate_rdpi(i$trait2, i$env),
    calculate_ESPI                = calculate_ESPI(i$trait2, i$env),
    calculate_espiid              = calculate_espiid(i$trait2),
    calculate_finlay_wilkinson    = calculate_finlay_wilkinson(i$gxe)
  ))
  for (nm in names(expected)) {
    expect_equal(got[[nm]], expected[[nm]], tolerance = 0, info = nm)
  }
})
