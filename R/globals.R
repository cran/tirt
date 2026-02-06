#' @import stats
#' @import utils
NULL

utils::globalVariables(c(
  "item", "model", "parameter", "response", "pvalue", "number",
  "Estimate", "Std_Error", "z_value", "Pr_z",
  "Estimate_slope", "Std_Error_slope", "z_value_slope", "Pr_z_slope",
  "Estimate_int", "parameter_int", "parameter_slope",
  "difficulty", "difficulty_se", "difficulty_zvalue", "difficulty_Pr_z",
  "discrimination", "discrimination_se", "discrimination_zvalue", "discrimination_Pr_z",
  "covariate", "covariate_se", "covariate_zvalue", "covariate_Pr_z",
  "infit", "outfit", "mean_theta", "mean_covariate",
  ".row_id", "theta", "step_1"
))
