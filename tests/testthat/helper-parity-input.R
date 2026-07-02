# Deterministic inputs shared by parity tests.
parity_inputs <- function() {
  set.seed(42)
  env <- 1:5
  trait <- c(2.0, 3.5, 5.0, 6.5, 8.0)          # clean linear
  trait2 <- c(1.0, 4.0, 2.0, 9.0, 3.0)          # noisy
  gxe <- matrix(c(2, 4, 6, 3, 6, 9, 1, 2, 3), nrow = 3, byrow = TRUE)  # 3 geno x 3 env
  list(env = env, trait = trait, trait2 = trait2, gxe = gxe)
}
