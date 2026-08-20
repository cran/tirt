#' Item Information Function
#'
#' @description
#' Computes the Fisher item information function for every item at a set of
#' ability (theta) points. The item information function shows how precisely an
#' item measures ability across the latent trait continuum: high information at
#' a given theta means the item discriminates well among examinees located
#' there. It is the building block of the test information function
#' (\code{\link{test_info}}) and of the conditional standard error of
#' measurement.
#'
#' @param item_params A data frame of item parameters, typically the
#'   \code{item_params} element returned by \code{\link{binary_irt}},
#'   \code{\link{polytomous_irt}}, or \code{\link{mixed_irt}}. The function
#'   auto-detects the model of each item from the available columns
#'   (discrimination, difficulty, guessing, thresholds/steps) or from a
#'   \code{model} column if present. All models in the package are supported:
#'   Rasch, 2PL, 3PL, GRM, GPCM, and PCM.
#' @param theta A numeric vector of ability values at which to evaluate the
#'   information (default = \code{seq(-4, 4, by = 0.1)}).
#' @param model Optional character string (applied to all items) or character
#'   vector of length equal to the number of items, forcing the model family
#'   used for each item. Advanced use only; when \code{NULL} (default) the model
#'   is detected automatically.
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}, matching the estimation functions in this package). Advanced
#'   users who calibrated on the normal-metric scale may set \code{D = 1.702}.
#'
#' @return A numeric matrix with one row per item and one column per theta point.
#'   Row names are the item names and column names are the theta values.
#'
#' @details
#' For dichotomous models the information is
#' \eqn{I(\theta) = (Da)^2 P(\theta)(1 - P(\theta))} for the Rasch and 2PL
#' models, and the standard Birnbaum form
#' \eqn{I(\theta) = (Da)^2 (1 - c)(1 - P^*)P^{*2} / P} for the 3PL model,
#' where \eqn{P^*} is the 2PL probability and \eqn{P} the full 3PL probability
#' (this reduces to the 2PL expression when the guessing parameter \eqn{c = 0}).
#' For polytomous models the general expected-information formula
#' \eqn{I(\theta) = \sum_k (\partial P_k / \partial \theta)^2 / P_k} is used,
#' which for the GPCM and PCM equals \eqn{(Da)^2 \mathrm{Var}(X \mid \theta)}.
#'
#' @seealso \code{\link{test_info}}, \code{\link{binary_irt}},
#'   \code{\link{polytomous_irt}}
#'
#' @examples
#'   # --- Example 1: dichotomous items ---
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 400,
#'                  item_structure = list(list(model = "2PL", n_items = 6)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   # Information for each item at a grid of abilities
#'   info <- item_info(fit$item_params, theta = seq(-3, 3, by = 0.5))
#'   round(info, 3)
#'
#'   \donttest{
#'   # --- Example 2: polytomous (GRM) items ---
#'   simp <- sim_irt(n_people = 400,
#'                   item_structure = list(list(model = "GRM", n_items = 5,
#'                                              categories = 4)))
#'   fitp <- polytomous_irt(simp$resp, model = "GRM", method = "EM",
#'                          control = list(max_iter = 15, verbose = FALSE))
#'   item_info(fitp$item_params, theta = c(-2, -1, 0, 1, 2))
#'   }
#' @export
item_info <- function(item_params,
                      theta = seq(-4, 4, by = 0.1),
                      model = NULL,
                      D = 1) {

  parsed <- .tirt_parse_params(item_params, model = model, D = D)
  theta <- as.numeric(theta)

  info_mat <- matrix(NA_real_, nrow = parsed$n_items, ncol = length(theta))
  for (j in seq_len(parsed$n_items)) {
    info_mat[j, ] <- .tirt_item_info(parsed$spec[[j]], theta)
  }

  rownames(info_mat) <- parsed$item
  colnames(info_mat) <- round(theta, 4)
  info_mat
}


#' Test Information Function and Conditional Standard Error
#'
#' @description
#' Computes the test information function (TIF) as the sum of all item
#' information values at each ability point, together with the conditional
#' standard error of measurement (SEM). The test information function
#' summarizes where along the ability scale the test measures most precisely,
#' and the conditional SEM (\eqn{1 / \sqrt{TIF}}) translates that precision back
#' onto the theta metric. This is a routine part of test evaluation and form
#' assembly in operational testing programs.
#'
#' @param item_params A data frame of item parameters (see \code{\link{item_info}}).
#' @param theta A numeric vector of ability values at which to evaluate the
#'   test information (default = \code{seq(-4, 4, by = 0.1)}).
#' @param model Optional model override (see \code{\link{item_info}}).
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}).
#'
#' @return A data frame with one row per theta point and the columns:
#' \itemize{
#'   \item \code{theta}: the ability value.
#'   \item \code{test_info}: the test information (sum of item information).
#'   \item \code{sem}: the conditional standard error of measurement,
#'     \eqn{1 / \sqrt{test\_info}}.
#'   \item \code{reliability}: the marginal-style conditional reliability
#'     \eqn{TIF / (TIF + 1)} for a standard-normal ability scale.
#' }
#'
#' @seealso \code{\link{item_info}}
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 400,
#'                  item_structure = list(list(model = "2PL", n_items = 10)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   tif <- test_info(fit$item_params, theta = seq(-3, 3, by = 0.5))
#'   print(tif)
#'
#'   # The ability where the test is most informative
#'   tif$theta[which.max(tif$test_info)]
#' @export
test_info <- function(item_params,
                      theta = seq(-4, 4, by = 0.1),
                      model = NULL,
                      D = 1) {

  theta <- as.numeric(theta)
  info_mat <- item_info(item_params, theta = theta, model = model, D = D)

  tif <- colSums(info_mat, na.rm = TRUE)
  sem <- ifelse(tif > 0, 1 / sqrt(tif), NA_real_)

  data.frame(
    theta       = theta,
    test_info   = tif,
    sem         = sem,
    reliability = tif / (tif + 1),
    row.names   = NULL,
    stringsAsFactors = FALSE
  )
}
