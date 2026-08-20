#' Person Fit Statistics (Standardized Log-Likelihood, lz)
#'
#' @description
#' Computes person-level fit statistics to detect examinees whose response
#' patterns are unlikely under the fitted item response model (for example,
#' careless responding, cheating, or unusual guessing). The function returns the
#' \eqn{l_z} standardized log-likelihood index of Drasgow, Levine, and Williams
#' (1985), which under model fit is approximately standard normal, and flags
#' persons whose \eqn{l_z} falls below a critical value.
#'
#' @param data A data frame of item responses (rows = persons, columns = items).
#'   Dichotomous items are 0/1; polytomous items are consecutive integers
#'   starting at 0. Missing values (\code{NA}) are allowed and are skipped
#'   person-by-person.
#' @param item_params A data frame of item parameters (the \code{item_params}
#'   element from \code{\link{binary_irt}}, \code{\link{polytomous_irt}}, or
#'   \code{\link{mixed_irt}}). Item order must match the columns of \code{data}.
#' @param theta A numeric vector of estimated person abilities (length equal to
#'   the number of rows in \code{data}), or a person-parameter data frame
#'   containing an \code{ability} or \code{theta} column (such as the
#'   \code{person_params} element returned by the estimation functions).
#' @param critical Numeric. Persons with \eqn{l_z} below this value are flagged
#'   as potentially misfitting (default = \code{-1.96}, the lower 2.5\% point of
#'   the standard normal distribution; \code{-2} is another common choice).
#' @param model Optional model override (see \code{\link{item_info}}).
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}).
#'
#' @return A data frame with one row per person and the columns:
#' \itemize{
#'   \item \code{person}: person index (row number).
#'   \item \code{n_items}: number of items answered.
#'   \item \code{theta}: the ability estimate used.
#'   \item \code{loglik}: the observed log-likelihood of the response pattern.
#'   \item \code{lz}: the standardized log-likelihood person-fit index.
#'   \item \code{flag}: logical; \code{TRUE} if \code{lz < critical}.
#' }
#'
#' @details
#' For a person with ability \eqn{\theta}, let \eqn{l_0} be the log-likelihood of
#' the observed responses. The index is
#' \eqn{l_z = (l_0 - E[l_0]) / \sqrt{\mathrm{Var}(l_0)}}, where the expectation
#' and variance are taken over the model-implied category probabilities at
#' \eqn{\theta}. Large negative values indicate response patterns that are less
#' likely than the model predicts.
#'
#' @references
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness
#' measurement with polychotomous item response models and standardized indices.
#' \emph{British Journal of Mathematical and Statistical Psychology, 38}(1),
#' 67-86.
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 300,
#'                  item_structure = list(list(model = "2PL", n_items = 15)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   pf <- person_fit(sim$resp, fit$item_params, fit$person_params)
#'   head(pf)
#'
#'   # How many examinees are flagged as misfitting?
#'   sum(pf$flag, na.rm = TRUE)
#' @export
person_fit <- function(data,
                       item_params,
                       theta,
                       critical = -1.96,
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

  th <- .tirt_extract_theta(theta)
  if (length(th) != N) {
    stop(sprintf("Length of 'theta' (%d) must match the number of rows in 'data' (%d).",
                 length(th), N))
  }

  # Evaluate item-by-item, reusing the per-item probability engine.
  # Probabilities are computed at the vector of person thetas once per item.
  obs_ll  <- numeric(N)
  exp_ll  <- numeric(N)
  var_ll  <- numeric(N)
  counts  <- integer(N)

  for (j in seq_len(J)) {
    sp <- parsed$spec[[j]]
    if (is.na(sp$a) || (sp$model %in% c("RASCH", "2PL", "3PL") && is.na(sp$b))) next

    probs <- .tirt_cat_probs(sp, th)          # N x n_cat
    probs <- pmin(pmax(probs, 1e-10), 1)
    logp  <- log(probs)

    x <- resp[, j]
    valid <- !is.na(x)
    if (!any(valid)) next

    # observed log-likelihood contribution
    idx <- cbind(which(valid), x[valid] + 1L)
    obs_ll[valid] <- obs_ll[valid] + logp[idx]

    # expectation and variance of the log-likelihood over categories
    e_term <- rowSums(probs * logp)                 # E[log P]
    v_term <- rowSums(probs * logp^2) - e_term^2    # Var[log P]
    exp_ll[valid] <- exp_ll[valid] + e_term[valid]
    var_ll[valid] <- var_ll[valid] + v_term[valid]
    counts[valid] <- counts[valid] + 1L
  }

  lz <- ifelse(var_ll > 0 & counts > 0,
               (obs_ll - exp_ll) / sqrt(var_ll),
               NA_real_)
  obs_ll[counts == 0] <- NA_real_

  out <- data.frame(
    person  = seq_len(N),
    n_items = counts,
    theta   = round(th, 3),
    loglik  = round(obs_ll, 3),
    lz      = round(lz, 3),
    flag    = lz < critical,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  out
}
