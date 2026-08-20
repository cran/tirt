#' Local Dependence Statistics (Yen's Q3)
#'
#' @description
#' Computes Yen's (1984) \eqn{Q_3} statistic for every item pair. \eqn{Q_3} is
#' the correlation between the model residuals of two items after the latent
#' trait has been partialled out. Under the local-independence assumption of
#' unidimensional item response theory these residual correlations should be
#' small and negative (around \eqn{-1/(J-1)}); large positive values signal that
#' two items share something beyond the common trait, such as a passage or a
#' scenario. \eqn{Q_3} is therefore used both to decide whether testlet modeling
#' is needed and to check whether a testlet model has adequately absorbed the
#' residual dependence.
#'
#' @param data A data frame of item responses (rows = persons, columns = items).
#'   Dichotomous items are 0/1; polytomous items are consecutive integers
#'   starting at 0.
#' @param item_params A data frame of item parameters (the \code{item_params}
#'   element from \code{\link{binary_irt}}, \code{\link{polytomous_irt}}, or
#'   \code{\link{mixed_irt}}). Item order must match the columns of \code{data}.
#' @param theta Optional numeric vector of estimated abilities (length equal to
#'   the number of rows in \code{data}), or a person-parameter data frame with
#'   an \code{ability}/\code{theta} column. If \code{NULL} (default), abilities
#'   are estimated internally by expected a posteriori (EAP) scoring from
#'   \code{item_params}.
#' @param theta_range Numeric vector of length 2 giving the ability grid bounds
#'   used for internal EAP scoring (default \code{c(-4, 4)}).
#' @param n_points Integer. Number of quadrature points for internal EAP scoring
#'   (default \code{41}).
#' @param model Optional model override (see \code{\link{item_info}}).
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}).
#'
#' @return A symmetric numeric matrix of \eqn{Q_3} statistics with one row and
#'   one column per item (diagonal set to \code{1}). The off-diagonal entries
#'   are the residual correlations for each item pair. The average off-diagonal
#'   \eqn{Q_3} is attached as the attribute \code{"mean_q3"} and the largest
#'   absolute value as \code{"max_abs_q3"}.
#'
#' @references
#' Yen, W. M. (1984). Effects of local item dependence on the fit and equating
#' performance of the three-parameter logistic model. \emph{Applied
#' Psychological Measurement, 8}(2), 125-145.
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 400,
#'                  item_structure = list(list(model = "2PL", n_items = 8)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   q3 <- ld_stats(sim$resp, fit$item_params)
#'   round(q3, 3)
#'   attr(q3, "max_abs_q3")
#' @export
ld_stats <- function(data,
                     item_params,
                     theta = NULL,
                     theta_range = c(-4, 4),
                     n_points = 41,
                     model = NULL,
                     D = 1) {

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix of responses.")
  }
  resp <- as.matrix(data)
  N <- nrow(resp)
  J <- ncol(resp)

  parsed <- .tirt_parse_params(item_params, model = model, D = D)
  if (parsed$n_items != J) {
    stop(sprintf("Number of items in 'item_params' (%d) does not match the number of columns in 'data' (%d).",
                 parsed$n_items, J))
  }
  specs <- parsed$spec

  # --- 1. Abilities: user-supplied or estimated by EAP ---
  if (is.null(theta)) {
    nodes <- seq(theta_range[1], theta_range[2], length.out = n_points)
    wts <- dnorm(nodes); wts <- wts / sum(wts)

    # log-likelihood of each person across nodes
    logL <- matrix(0, N, n_points)
    for (j in seq_len(J)) {
      sp <- specs[[j]]
      if (is.na(sp$a) || (sp$model %in% c("RASCH", "2PL", "3PL") && is.na(sp$b))) next
      pr <- pmin(pmax(.tirt_cat_probs(sp, nodes), 1e-10), 1)   # nodes x n_cat
      x <- resp[, j]
      valid <- !is.na(x)
      if (!any(valid)) next
      lp <- log(pr)                                            # nodes x n_cat
      logL[valid, ] <- logL[valid, ] + t(lp[, x[valid] + 1L, drop = FALSE])
    }
    Fpost <- exp(logL - apply(logL, 1, max)) * matrix(wts, N, n_points, byrow = TRUE)
    Fpost <- Fpost / rowSums(Fpost)
    th <- as.vector(Fpost %*% nodes)
  } else {
    th <- .tirt_extract_theta(theta)
    if (length(th) != N) {
      stop(sprintf("Length of 'theta' (%d) must match the number of rows in 'data' (%d).",
                   length(th), N))
    }
  }

  # --- 2. Residuals: observed minus expected item score ---
  resid <- matrix(NA_real_, N, J)
  for (j in seq_len(J)) {
    sp <- specs[[j]]
    if (is.na(sp$a) || (sp$model %in% c("RASCH", "2PL", "3PL") && is.na(sp$b))) next
    expct <- .tirt_expected_score(sp, th)
    resid[, j] <- resp[, j] - expct
  }

  # --- 3. Pairwise correlations of residuals ---
  q3 <- suppressWarnings(cor(resid, use = "pairwise.complete.obs"))
  q3 <- as.matrix(q3)
  rownames(q3) <- colnames(q3) <- parsed$item
  diag(q3) <- 1

  off <- q3[upper.tri(q3)]
  off <- off[is.finite(off)]
  attr(q3, "mean_q3") <- if (length(off)) mean(off) else NA_real_
  attr(q3, "max_abs_q3") <- if (length(off)) max(abs(off)) else NA_real_
  q3
}
