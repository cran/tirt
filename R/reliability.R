#' Test Reliability Indices
#'
#' @description
#' Computes commonly reported reliability coefficients for a test. Depending on
#' the information supplied it returns the empirical (marginal) reliability from
#' latent-trait estimates and their standard errors, the model-based marginal
#' reliability from the test information function, and the classical Cronbach's
#' alpha from the raw responses. Reliability summarizes how consistently the test
#' orders examinees and is a standard entry in every technical report.
#'
#' @param person_params Optional. A person-parameter data frame containing an
#'   \code{ability}/\code{theta} column and a standard-error column (such as
#'   \code{ability_se} or \code{se}); typically the \code{person_params} element
#'   returned by the estimation functions. Used for the empirical reliability.
#' @param data Optional. A data frame of item responses used to compute
#'   Cronbach's alpha (rows = persons, columns = items).
#' @param item_params Optional. A data frame of item parameters used to compute
#'   the model-based marginal reliability from the test information function.
#' @param prior_sd Numeric. Standard deviation of the population ability
#'   distribution used for the model-based marginal reliability (default
#'   \code{1}).
#' @param theta_range Numeric vector of length 2 giving the ability grid bounds
#'   for the marginal-reliability integration (default \code{c(-4, 4)}).
#' @param n_points Integer. Number of grid points for the marginal-reliability
#'   integration (default \code{61}).
#' @param model Optional model override (see \code{\link{item_info}}).
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}).
#'
#' @return A data frame with the columns \code{Index} and \code{Value},
#'   reporting whichever reliability coefficients could be computed from the
#'   supplied arguments:
#' \itemize{
#'   \item \code{empirical_reliability}: from person estimates and their SEs.
#'   \item \code{marginal_reliability}: from the test information function.
#'   \item \code{cronbach_alpha}: classical internal-consistency reliability.
#'   \item \code{n_items}, \code{n_persons}: sample descriptors when available.
#' }
#'
#' @details
#' The empirical reliability is
#' \eqn{\mathrm{Var}(\hat\theta) / (\mathrm{Var}(\hat\theta) + \overline{SE^2})}.
#' The model-based marginal reliability is
#' \eqn{\sigma^2 / (\sigma^2 + \overline{1/I(\theta)})}, where the average error
#' variance is taken over the population ability distribution. Cronbach's alpha
#' is \eqn{\frac{J}{J-1}\left(1 - \frac{\sum_j s_j^2}{s_T^2}\right)}.
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 500,
#'                  item_structure = list(list(model = "2PL", n_items = 12)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   reliability(person_params = fit$person_params,
#'               data = sim$resp,
#'               item_params = fit$item_params)
#' @export
reliability <- function(person_params = NULL,
                        data = NULL,
                        item_params = NULL,
                        prior_sd = 1,
                        theta_range = c(-4, 4),
                        n_points = 61,
                        model = NULL,
                        D = 1) {

  idx <- character(0)
  val <- numeric(0)
  add <- function(name, value) { idx[[length(idx) + 1]] <<- name; val[[length(val) + 1]] <<- value }

  # --- Empirical reliability from person estimates ---
  if (!is.null(person_params)) {
    th <- .tirt_extract_theta(person_params)
    se <- NULL
    if (is.data.frame(person_params)) {
      for (cand in c("ability_se", "se", "theta_se", "SE", "SE_Theta")) {
        if (cand %in% names(person_params)) { se <- as.numeric(person_params[[cand]]); break }
      }
    }
    if (!is.null(se)) {
      ok <- is.finite(th) & is.finite(se)
      if (sum(ok) > 1) {
        v_th <- var(th[ok])
        m_se2 <- mean(se[ok]^2)
        add("empirical_reliability", round(v_th / (v_th + m_se2), 4))
      }
    }
  }

  # --- Marginal reliability from test information ---
  if (!is.null(item_params)) {
    nodes <- seq(theta_range[1], theta_range[2], length.out = n_points)
    w <- dnorm(nodes, 0, prior_sd); w <- w / sum(w)
    tif <- test_info(item_params, theta = nodes, model = model, D = D)$test_info
    err_var <- sum(w * ifelse(tif > 0, 1 / tif, NA_real_), na.rm = TRUE)
    sig2 <- prior_sd^2
    add("marginal_reliability", round(sig2 / (sig2 + err_var), 4))
  }

  # --- Classical Cronbach's alpha ---
  if (!is.null(data)) {
    mat <- as.matrix(data)
    cc <- mat[stats::complete.cases(mat), , drop = FALSE]
    J <- ncol(cc)
    if (nrow(cc) > 1 && J > 1) {
      item_var <- apply(cc, 2, var)
      total_var <- var(rowSums(cc))
      if (total_var > 0) {
        alpha <- (J / (J - 1)) * (1 - sum(item_var) / total_var)
        add("cronbach_alpha", round(alpha, 4))
      }
      add("n_items", J)
      add("n_persons", nrow(cc))
    }
  }

  if (length(idx) == 0) {
    stop("Nothing to compute: supply at least one of 'person_params', 'item_params', or 'data'.")
  }

  data.frame(Index = unlist(idx), Value = unlist(val),
             row.names = NULL, stringsAsFactors = FALSE)
}
