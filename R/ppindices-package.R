#' ppindices: Phenotypic Plasticity Indices
#'
#' A collection of implementations for quantifying phenotypic plasticity from
#' reaction-norm data. Provides functions to compute a broad range of
#' plasticity indices from trait values measured across environments.
#'
#' @keywords internal
#' @importFrom stats AIC BIC aggregate anova aov as.formula coef cor cov dist
#'   lm median prcomp predict quantile sd setNames var
#' @importFrom graphics abline legend points
#' @importFrom utils combn install.packages setTxtProgressBar txtProgressBar
"_PACKAGE"

## Names that appear as bare (unquoted) symbols in non-standard evaluation
## contexts, or that are referenced without being defined (pre-existing bugs
## in the ported reference implementations; see calculate_CVm() and
## calculate_reaction_norm_non_linear() documentation). Declared here only to
## silence R CMD check's "no visible binding for global variable" NOTE; the
## underlying behaviour of these functions is intentionally left unchanged.
utils::globalVariables(c("env", "group_labels", "Combined_Factors", "Trait_Value"))
