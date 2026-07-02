test_that("every reference index is exported by ppindices", {
  expected <- c(
    "calculate_CVt","calculate_reaction_norm_slope","calculate_reaction_norm_non_linear",
    "calculate_D_slope","calculate_RC","calculate_CVm","calculate_CVmd",
    "calculate_grand_plasticity","calculate_PPF","calculate_Phenotypic_Plasticity_Index",
    "calculate_PImd","calculate_PILSM","calculate_RTR","calculate_PIR",
    "calculate_rdpi","rdpi_mean_calculation","calculate_ESPI","calculate_espiid",
    "evwpi_calculation","calculate_PSI","calculate_RPI","calculate_PQ","calculate_PR",
    "calculate_NRW","calculate_ESP","calculate_general_PD","calculate_FPI","calculate_TPS",
    "calculate_DPI","calculate_CEV","calculate_PRI","calculate_PFI","calculate_SPI",
    "calculate_APC","calculate_SI","calculate_RSI","calculate_EVS","calculate_MVPi",
    "calculate_SPM","calculate_Plasticity_Ratio","calculate_plasticity","cross_env_cov",
    "calculate_finlay_wilkinson"
  )
  exported <- getNamespaceExports("ppindices")
  missing <- setdiff(expected, exported)
  expect_identical(missing, character(0))
})
