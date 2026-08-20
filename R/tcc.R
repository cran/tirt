#' Test and Item Characteristic Curves (Expected Scores)
#'
#' @description
#' Computes item characteristic curves (the expected score of each item as a
#' function of ability) and the test characteristic curve (the expected total
#' score as a function of ability). The test characteristic curve maps the latent
#' ability scale onto the number-correct (true-score) scale and underlies
#' true-score equating and score reporting; the item curves show how each item's
#' expected score rises across the trait.
#'
#' @param item_params A data frame of item parameters (the \code{item_params}
#'   element from \code{\link{binary_irt}}, \code{\link{polytomous_irt}}, or
#'   \code{\link{mixed_irt}}). All models in the package are supported.
#' @param theta A numeric vector of ability values at which to evaluate the
#'   curves (default = \code{seq(-4, 4, by = 0.1)}).
#' @param model Optional model override (see \code{\link{item_info}}).
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}).
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{test_curve}: a data frame with columns \code{theta} and
#'     \code{expected_score} (the test characteristic curve).
#'   \item \code{item_curves}: a numeric matrix with one row per item and one
#'     column per theta value, giving the expected score of each item.
#' }
#'
#' @seealso \code{\link{item_info}}, \code{\link{score_table}}
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 400,
#'                  item_structure = list(list(model = "2PL", n_items = 10)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   curves <- tcc(fit$item_params, theta = seq(-3, 3, by = 0.5))
#'   curves$test_curve
#'
#'   # Expected score on the first item across ability
#'   curves$item_curves[1, ]
#' @export
tcc <- function(item_params,
                theta = seq(-4, 4, by = 0.1),
                model = NULL,
                D = 1) {

  parsed <- .tirt_parse_params(item_params, model = model, D = D)
  theta <- as.numeric(theta)

  es_mat <- matrix(NA_real_, nrow = parsed$n_items, ncol = length(theta))
  for (j in seq_len(parsed$n_items)) {
    sp <- parsed$spec[[j]]
    if (is.na(sp$a) || (sp$model %in% c("RASCH", "2PL", "3PL") && is.na(sp$b))) next
    es_mat[j, ] <- .tirt_expected_score(sp, theta)
  }
  rownames(es_mat) <- parsed$item
  colnames(es_mat) <- round(theta, 4)

  test_curve <- data.frame(
    theta = theta,
    expected_score = colSums(es_mat, na.rm = TRUE),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  list(test_curve = test_curve, item_curves = es_mat)
}
