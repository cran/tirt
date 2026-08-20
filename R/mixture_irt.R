#' Mixture Item Response Theory Model (Latent-Class Rasch / 2PL)
#'
#' @description
#' Fits a mixture item response theory model (Rost, 1990) in which the examinee
#' population is assumed to consist of a small number of unobserved latent
#' classes, each with its own set of item parameters. Mixture IRT models are used
#' to detect qualitatively different response strategies, unmodeled
#' subpopulations, or classes for which item difficulty ordering differs. The
#' model is estimated by the Expectation-Maximization algorithm with a fixed
#' standard-normal ability distribution within each class.
#'
#' @param data A data frame of dichotomous (0/1) item responses (rows = persons,
#'   columns = items).
#' @param n_class Integer. Number of latent classes to estimate (default
#'   \code{2}).
#' @param model String. \code{"Rasch"} (default) or \code{"2PL"}; controls
#'   whether class-specific discriminations are estimated.
#' @param control A \code{list} of control parameters for the algorithm:
#'   \itemize{
#'     \item \code{max_iter}: Maximum number of EM iterations (default = 100).
#'     \item \code{converge_tol}: Convergence criterion on the log-likelihood
#'       change (default = 1e-4).
#'     \item \code{quad_points}: Number of quadrature points for the within-class
#'       ability distribution (default = 21).
#'     \item \code{theta_range}: Ability grid bounds (default = c(-4, 4)).
#'     \item \code{verbose}: Logical; print progress (default = TRUE).
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{item_params}: A data frame of class-specific item difficulties
#'     (and discriminations for the 2PL model), with one column per class.
#'   \item \code{class_params}: A data frame of estimated mixing proportions for
#'     each latent class.
#'   \item \code{person_params}: A data frame of posterior class-membership
#'     probabilities and the modal (most likely) class for each person.
#'   \item \code{model_fit}: A data frame with the log-likelihood, AIC, BIC,
#'     number of classes, and classification entropy.
#' }
#'
#' @references
#' Rost, J. (1990). Rasch models in latent classes: An integration of two
#' approaches to item analysis. \emph{Applied Psychological Measurement, 14}(3),
#' 271-282.
#'
#' @examples
#'   # Two classes with reversed difficulty ordering
#'   set.seed(2025)
#'   N <- 300; J <- 8
#'   b1 <- seq(-1.5, 1.5, length.out = J)
#'   b2 <- rev(b1)
#'   theta <- rnorm(N)
#'   cls <- rep(1:2, each = N / 2)
#'   resp <- matrix(0, N, J)
#'   for (i in 1:N) {
#'     b <- if (cls[i] == 1) b1 else b2
#'     resp[i, ] <- rbinom(J, 1, 1 / (1 + exp(-(theta[i] - b))))
#'   }
#'   df <- as.data.frame(resp); names(df) <- paste0("I", 1:J)
#'
#'   fit <- mixture_irt(df, n_class = 2, model = "Rasch",
#'                      control = list(max_iter = 30, verbose = FALSE))
#'   fit$class_params
#'   head(fit$item_params)
#' @export
mixture_irt <- function(data,
                        n_class = 2,
                        model = "Rasch",
                        control = list()) {

  con <- list(max_iter = 100, converge_tol = 1e-4, quad_points = 21,
              theta_range = c(-4, 4), verbose = TRUE)
  con[names(control)] <- control

  model <- toupper(model)
  if (!model %in% c("RASCH", "2PL")) stop("'model' must be 'Rasch' or '2PL'.")
  if (!is.data.frame(data) && !is.matrix(data)) stop("'data' must be a data frame or matrix.")

  X <- as.matrix(data)
  if (!all(X[!is.na(X)] %in% c(0, 1))) stop("mixture_irt() currently supports dichotomous (0/1) data.")
  N <- nrow(X); J <- ncol(X)
  C <- as.integer(n_class)
  item_names <- colnames(X); if (is.null(item_names)) item_names <- paste0("Item_", seq_len(J))

  nodes <- seq(con$theta_range[1], con$theta_range[2], length.out = con$quad_points)
  qw <- dnorm(nodes); qw <- qw / sum(qw)
  Q <- length(nodes)

  # --- initialize ---
  pi_c <- rep(1 / C, C)
  p_mean <- pmin(pmax(colMeans(X, na.rm = TRUE), 0.05), 0.95)
  b0 <- -log(p_mean / (1 - p_mean))
  # spread classes apart with small perturbations for identifiability
  b <- matrix(rep(b0, C), J, C)
  for (cc in seq_len(C)) b[, cc] <- b0 + seq(-0.5, 0.5, length.out = J) * (cc - (C + 1) / 2)
  a <- matrix(1, J, C)

  prob_node <- function(theta, aj, bj) {
    z <- aj * (theta - bj); z <- pmin(pmax(z, -30), 30); 1 / (1 + exp(-z))
  }

  if (con$verbose) {
    cat(sprintf("\nStarting Mixture %s Estimation (%d classes)...\n", model, C))
    cat("------------------------------------------------\n")
  }

  old_ll <- -Inf
  is_converged <- FALSE
  Post_class <- matrix(1 / C, N, C)

  for (iter in seq_len(con$max_iter)) {

    # E-step: class x node likelihoods
    # L_node[i, q, c] then marginal over q
    Lc <- matrix(0, N, C)                    # marginal per class
    node_post <- array(0, dim = c(N, Q, C))  # unnormalized person-node-class

    for (cc in seq_len(C)) {
      logLq <- matrix(0, N, Q)
      for (j in seq_len(J)) {
        P <- prob_node(nodes, a[j, cc], b[j, cc])
        P <- pmin(pmax(P, 1e-10), 1 - 1e-10)
        x <- X[, j]; valid <- !is.na(x)
        if (!any(valid)) next
        add <- outer(x[valid], log(P)) + outer(1 - x[valid], log(1 - P))
        logLq[valid, ] <- logLq[valid, ] + add
      }
      m <- apply(logLq, 1, max)
      Fq <- exp(logLq - m) * matrix(qw, N, Q, byrow = TRUE)
      Lc[, cc] <- rowSums(Fq) * exp(m)
      node_post[, , cc] <- Fq * exp(m)       # unnormalized joint of node & data
    }

    joint <- sweep(Lc, 2, pi_c, "*")         # N x C, pi_c * L_ic
    denom <- rowSums(joint)
    denom[denom <= 0] <- 1e-300
    Post_class <- joint / denom
    ll <- sum(log(denom))

    # M-step: mixing proportions
    pi_c <- colMeans(Post_class)
    pi_c <- pmax(pi_c, 1e-6); pi_c <- pi_c / sum(pi_c)

    # M-step: item parameters per class
    for (cc in seq_len(C)) {
      # person-node weights within class cc: r_ic * node_post / Lc
      denom_c <- Lc[, cc]; denom_c[denom_c <= 0] <- 1e-300
      w_node <- (node_post[, , cc] / denom_c) * Post_class[, cc]   # N x Q

      for (j in seq_len(J)) {
        x <- X[, j]; valid <- !is.na(x)
        wj <- w_node[valid, , drop = FALSE]
        r_q <- colSums(wj * x[valid])
        n_q <- colSums(wj)
        if (sum(n_q) < 1e-6) next

        # Weighted (expected-count) log-likelihood for this item in this class.
        ll_item <- function(par) {
          la <- if (model == "2PL") par[1] else 1
          lb <- if (model == "2PL") par[2] else par[1]
          P <- prob_node(nodes, la, lb); P <- pmin(pmax(P, 1e-9), 1 - 1e-9)
          sum(r_q * log(P) + (n_q - r_q) * log(1 - P))
        }

        # Damped Newton-Raphson with finite-difference gradient / diagonal Hessian.
        par <- if (model == "2PL") c(a[j, cc], b[j, cc]) else b[j, cc]
        h <- 1e-4
        for (nr in 1:15) {
          f0 <- ll_item(par)
          grad <- numeric(length(par)); hess <- numeric(length(par))
          for (k in seq_along(par)) {
            pu <- par; pu[k] <- pu[k] + h
            pd <- par; pd[k] <- pd[k] - h
            grad[k] <- (ll_item(pu) - ll_item(pd)) / (2 * h)
            hess[k] <- (ll_item(pu) - 2 * f0 + ll_item(pd)) / (h^2)
          }
          hess <- ifelse(hess > -1e-6, -1e-6, hess)   # keep negative definite
          step <- grad / hess
          step <- pmax(pmin(step, 1), -1)
          par <- par - step
          if (model == "2PL") par[1] <- max(0.05, min(4, par[1]))
          par[length(par)] <- max(-10, min(10, par[length(par)]))
          if (max(abs(step)) < 1e-4) break
        }
        if (model == "2PL") { a[j, cc] <- par[1]; b[j, cc] <- par[2] } else b[j, cc] <- par[1]
      }
    }

    if (con$verbose) cat(sprintf("\rIteration %d: LogLik = %.3f   ", iter, ll))

    if (iter > 2 && abs(ll - old_ll) < con$converge_tol) {
      if (con$verbose) cat("\n\n>>> Convergence Confirmation: Model Converged!\n")
      is_converged <- TRUE
      break
    }
    old_ll <- ll
  }

  if (!is_converged && con$verbose) {
    cat(sprintf("\n\n>>> NOTE: Stopped because max Iteration Reached at %d.\n", con$max_iter))
  }

  # --- Output ---
  item_df <- data.frame(item = item_names, stringsAsFactors = FALSE)
  for (cc in seq_len(C)) {
    if (model == "2PL") item_df[[paste0("discrimination_c", cc)]] <- round(a[, cc], 3)
    item_df[[paste0("difficulty_c", cc)]] <- round(b[, cc], 3)
  }

  class_df <- data.frame(class = seq_len(C), proportion = round(pi_c, 4),
                         stringsAsFactors = FALSE)

  modal <- max.col(Post_class, ties.method = "first")
  person_df <- data.frame(person = seq_len(N))
  for (cc in seq_len(C)) person_df[[paste0("class", cc)]] <- round(Post_class[, cc], 3)
  person_df$modal_class <- modal

  # classification entropy (normalized); undefined for a single class
  ent <- -sum(Post_class * log(pmax(Post_class, 1e-12)))
  entropy <- if (C <= 1) NA_real_ else 1 - ent / (N * log(C))

  n_par <- (C - 1) + C * J * (if (model == "2PL") 2 else 1)
  fit_df <- data.frame(
    Index = c("LogLikelihood", "AIC", "BIC", "n_class", "entropy"),
    Value = c(round(old_ll, 3),
              round(2 * n_par - 2 * old_ll, 3),
              round(n_par * log(N) - 2 * old_ll, 3),
              C, round(entropy, 3)),
    stringsAsFactors = FALSE
  )

  if (con$verbose) {
    cat("Finished All Estimation.\n")
    cat("------------------------------------------------\n")
  }

  list(item_params = item_df,
       class_params = class_df,
       person_params = person_df,
       model_fit = fit_df)
}
