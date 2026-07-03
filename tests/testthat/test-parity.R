test_that("packaged indices reproduce reference outputs bit-for-bit", {
  # Golden values captured from the original implementations in
  # CRC_1644_Z2_GWAS_simple/R-files/Plasticity_scores/ (see helper for the
  # shared inputs and call map). Every exported calculate_* index is covered.
  expected <- readRDS(test_path("fixtures", "parity.rds"))
  calls <- parity_calls(getNamespace("ppindices"))

  # The fixture defines the comparable set; every entry must have a call.
  missing <- setdiff(names(expected), names(calls))
  expect_identical(missing, character(0),
                   info = paste("no call for:", paste(missing, collapse = ", ")))

  for (nm in names(expected)) {
    got <- suppressWarnings(suppressMessages(calls[[nm]]()))
    expect_equal(got, expected[[nm]], tolerance = 0, info = nm)
  }
})
