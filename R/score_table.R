# Internal: Lord-Wingersky recursion.
# Given a list of item category-probability matrices (each nodes x n_cat),
# returns a matrix (nodes x (max_total_score + 1)) of the probability of every
# possible summed score at each ability node.
.tirt_lw_recursion <- function(prob_list) {
  n_nodes <- nrow(prob_list[[1]])
  # start: one item with a single "score 0" mass of 1
  dist <- matrix(1, n_nodes, 1)          # P(total = 0) = 1 before any item
  for (p in prob_list) {
    n_cat <- ncol(p)
    cur_max <- ncol(dist) - 1
    new_max <- cur_max + (n_cat - 1)
    new_dist <- matrix(0, n_nodes, new_max + 1)
    for (k in seq_len(n_cat)) {          # category score = k - 1
      s <- k - 1
      new_dist[, (1 + s):(cur_max + 1 + s)] <-
        new_dist[, (1 + s):(cur_max + 1 + s)] + dist * p[, k]
    }
    dist <- new_dist
  }
  dist
}

#' Summed-Score to Theta Conversion Table
#'
#' @description
#' Generates a raw-score (summed-score) to scale-score conversion table, mapping
#' every possible integer summed score to its corresponding ability estimate and
#' standard error. Such tables are standard in operational testing programs
#' because they let each examinee's number-correct (or number-of-points) score
#' be translated directly into a theta estimate without needing the full
#' response pattern. Three scoring methods are supported: expected a posteriori
#' (\code{"EAP"}), weighted likelihood (\code{"WLE"}), and maximum likelihood
#' (\code{"MLE"}).
#'
#' @param item_params A data frame of item parameters (the \code{item_params}
#'   element from \code{\link{binary_irt}}, \code{\link{polytomous_irt}}, or
#'   \code{\link{mixed_irt}}). All models in the package are supported.
#' @param method Character. Scoring method: \code{"EAP"} (default),
#'   \code{"WLE"}, or \code{"MLE"}.
#' @param prior_mean,prior_sd Mean and standard deviation of the normal
#'   population prior used for the EAP method and for the quadrature weighting
#'   (defaults \code{0} and \code{1}).
#' @param theta_range Numeric vector of length 2 giving the ability grid bounds
#'   (default \code{c(-4, 4)}).
#' @param n_points Integer. Number of ability grid points used to evaluate the
#'   summed-score distributions (default \code{81}).
#' @param model Optional model override (see \code{\link{item_info}}).
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}).
#'
#' @return A data frame with one row per possible summed score and the columns:
#' \itemize{
#'   \item \code{summed_score}: the integer summed score (0 to the maximum
#'     possible score).
#'   \item \code{theta}: the ability estimate for that summed score.
#'   \item \code{se}: the standard error of the ability estimate.
#' }
#'
#' @details
#' The distribution of the summed score at each ability grid point is obtained
#' with the Lord-Wingersky (1984) recursion, which works for any mixture of
#' dichotomous and polytomous items. The \code{"EAP"} estimate is the mean of
#' the posterior of theta given the summed score, with \code{se} equal to the
#' posterior standard deviation. The \code{"MLE"} estimate maximizes the
#' summed-score likelihood; the \code{"WLE"} estimate (Warm, 1989) maximizes the
#' summed-score likelihood weighted by the square root of the test information,
#' which reduces the bias of the MLE. At the two extreme summed scores the MLE is
#' not finite and is reported at the boundary of \code{theta_range} with a missing
#' standard error; the WLE remains finite there and is estimated normally.
#'
#' @references
#' Lord, F. M., & Wingersky, M. S. (1984). Comparison of IRT true-score and
#' equipercentile observed-score equatings. \emph{Applied Psychological
#' Measurement, 8}(4), 453-461.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. \emph{Psychometrika, 54}(3), 427-450.
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 400,
#'                  item_structure = list(list(model = "2PL", n_items = 10)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   # EAP conversion table (0 to 10 correct)
#'   score_table(fit$item_params, method = "EAP")
#'
#'   \donttest{
#'   # Maximum-likelihood conversion table
#'   score_table(fit$item_params, method = "MLE")
#'   }
#' @export
score_table <- function(item_params,
                        method = c("EAP", "WLE", "MLE"),
                        prior_mean = 0,
                        prior_sd = 1,
                        theta_range = c(-4, 4),
                        n_points = 81,
                        model = NULL,
                        D = 1) {

  method <- match.arg(method)
  parsed <- .tirt_parse_params(item_params, model = model, D = D)

  # Drop items whose parameters are missing (0-variance items etc.)
  keep <- vapply(parsed$spec, function(s) {
    !is.na(s$a) && !(s$model %in% c("RASCH", "2PL", "3PL") && is.na(s$b))
  }, logical(1))
  specs <- parsed$spec[keep]
  if (length(specs) == 0) stop("No usable items found in 'item_params'.")

  max_score <- sum(vapply(specs, .tirt_max_score, integer(1)))

  nodes <- seq(theta_range[1], theta_range[2], length.out = n_points)
  prior <- dnorm(nodes, prior_mean, prior_sd)
  prior <- prior / sum(prior)

  # Summed-score distributions at each node (nodes x (max_score + 1))
  prob_list <- lapply(specs, .tirt_cat_probs, theta = nodes)
  lw <- .tirt_lw_recursion(prob_list)

  # Test information at each node (for the WLE weight)
  info_mat <- vapply(specs, .tirt_item_info, numeric(length(nodes)), theta = nodes)
  tif <- rowSums(matrix(info_mat, nrow = length(nodes)))

  scores <- 0:max_score
  theta_hat <- rep(NA_real_, length(scores))
  se_hat    <- rep(NA_real_, length(scores))

  for (r in scores) {
    lik <- lw[, r + 1]                                  # P(sum = r | node)

    if (method == "EAP") {
      post <- lik * prior
      denom <- sum(post)
      if (denom <= 0) next
      post <- post / denom
      m  <- sum(post * nodes)
      v  <- sum(post * (nodes - m)^2)
      theta_hat[r + 1] <- m
      se_hat[r + 1]    <- sqrt(v)

    } else {
      # MLE / WLE: maximize a (weighted) log-likelihood over the grid
      logL <- log(pmax(lik, 1e-300))
      if (method == "WLE") logL <- logL + 0.5 * log(pmax(tif, 1e-300))

      # The MLE diverges at the two extreme scores (monotone likelihood); Warm's
      # WLE stays finite there, so only the MLE is pinned to the boundary.
      if (method == "MLE" && (r == 0 || r == max_score)) {
        theta_hat[r + 1] <- if (r == 0) theta_range[1] else theta_range[2]
        se_hat[r + 1]    <- NA_real_
        next
      }

      best <- which.max(logL)
      # parabolic refinement around the grid maximum
      if (best > 1 && best < length(nodes)) {
        y1 <- logL[best - 1]; y2 <- logL[best]; y3 <- logL[best + 1]
        denom <- (y1 - 2 * y2 + y3)
        h <- nodes[2] - nodes[1]
        shift <- if (denom != 0) 0.5 * (y1 - y3) / denom else 0
        shift <- max(min(shift, 1), -1)
        theta_hat[r + 1] <- nodes[best] + shift * h
        # curvature -> standard error
        curv <- -denom / (h^2)
        se_hat[r + 1] <- if (curv > 0) 1 / sqrt(curv) else NA_real_
      } else {
        theta_hat[r + 1] <- nodes[best]
        se_hat[r + 1]    <- NA_real_
      }
    }
  }

  data.frame(
    summed_score = scores,
    theta        = round(theta_hat, 3),
    se           = round(se_hat, 3),
    row.names    = NULL,
    stringsAsFactors = FALSE
  )
}
