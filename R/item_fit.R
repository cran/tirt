#' Item Fit Statistics (Infit and Outfit)
#'
#' @description
#' Computes item-level fit statistics based on standardized response residuals:
#' the information-weighted mean square (infit) and the unweighted mean square
#' (outfit), together with their standardized (ZSTD) transformations. These
#' statistics flag items whose observed responses depart from the fitted model,
#' either because of noise in the tails of the ability distribution (outfit) or
#' near the item location (infit). The statistics are computed for dichotomous
#' and polytomous items alike.
#'
#' @param data A data frame of item responses (rows = persons, columns = items).
#'   Dichotomous items are 0/1; polytomous items are consecutive integers
#'   starting at 0. Item order must match \code{item_params}.
#' @param item_params A data frame of item parameters (the \code{item_params}
#'   element from \code{\link{binary_irt}}, \code{\link{polytomous_irt}}, or
#'   \code{\link{mixed_irt}}).
#' @param theta Optional numeric vector of estimated abilities, or a person-
#'   parameter data frame with an \code{ability}/\code{theta} column. If
#'   \code{NULL} (default), abilities are estimated internally by EAP scoring.
#' @param theta_range Numeric vector of length 2 giving the ability grid bounds
#'   used for internal EAP scoring (default \code{c(-4, 4)}).
#' @param n_points Integer. Number of quadrature points for internal EAP scoring
#'   (default \code{41}).
#' @param model Optional model override (see \code{\link{item_info}}).
#' @param D Scaling constant for the dichotomous logistic metric (default
#'   \code{1}).
#'
#' @return A data frame with one row per item and the columns:
#' \itemize{
#'   \item \code{item}: item name.
#'   \item \code{n}: number of responses used.
#'   \item \code{outfit}: outfit mean square (unweighted).
#'   \item \code{outfit_z}: standardized outfit (ZSTD).
#'   \item \code{infit}: infit mean square (information-weighted).
#'   \item \code{infit_z}: standardized infit (ZSTD).
#' }
#' Mean-square values near 1 indicate good fit; values above about 1.3-1.5
#' suggest underfit (noise), and values below about 0.7 suggest overfit
#' (dependency or redundancy).
#'
#' @references
#' Wright, B. D., & Masters, G. N. (1982). \emph{Rating scale analysis}. MESA
#' Press.
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 400,
#'                  item_structure = list(list(model = "2PL", n_items = 8)))
#'   fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
#'                     control = list(max_iter = 15, verbose = FALSE))
#'
#'   item_fit(sim$resp, fit$item_params)
#' @export
item_fit <- function(data,
                     item_params,
                     theta = NULL,
                     theta_range = c(-4, 4),
                     n_points = 41,
                     model = NULL,
                     D = 1) {

  if (!is.data.frame(data) && !is.matrix(data)) stop("'data' must be a data frame or matrix.")
  resp <- as.matrix(data)
  N <- nrow(resp); J <- ncol(resp)

  parsed <- .tirt_parse_params(item_params, model = model, D = D)
  if (parsed$n_items != J) {
    stop(sprintf("Number of items in 'item_params' (%d) does not match the number of columns in 'data' (%d).",
                 parsed$n_items, J))
  }
  specs <- parsed$spec

  # --- abilities: supplied or EAP ---
  if (is.null(theta)) {
    nodes <- seq(theta_range[1], theta_range[2], length.out = n_points)
    wts <- dnorm(nodes); wts <- wts / sum(wts)
    logL <- matrix(0, N, n_points)
    for (j in seq_len(J)) {
      sp <- specs[[j]]
      if (is.na(sp$a) || (sp$model %in% c("RASCH", "2PL", "3PL") && is.na(sp$b))) next
      pr <- pmin(pmax(.tirt_cat_probs(sp, nodes), 1e-10), 1)
      x <- resp[, j]; valid <- !is.na(x)
      if (!any(valid)) next
      logL[valid, ] <- logL[valid, ] + t(log(pr)[, x[valid] + 1L, drop = FALSE])
    }
    Fpost <- exp(logL - apply(logL, 1, max)) * matrix(wts, N, n_points, byrow = TRUE)
    Fpost <- Fpost / rowSums(Fpost)
    th <- as.vector(Fpost %*% nodes)
  } else {
    th <- .tirt_extract_theta(theta)
    if (length(th) != N) stop("Length of 'theta' must match the number of rows in 'data'.")
  }

  out <- vector("list", J)
  for (j in seq_len(J)) {
    sp <- specs[[j]]
    x <- resp[, j]; valid <- !is.na(x)
    if (is.na(sp$a) || (sp$model %in% c("RASCH", "2PL", "3PL") && is.na(sp$b)) || !any(valid)) {
      out[[j]] <- data.frame(item = parsed$item[j], n = sum(valid),
                             outfit = NA, outfit_z = NA, infit = NA, infit_z = NA,
                             stringsAsFactors = FALSE)
      next
    }

    thj <- th[valid]; xj <- x[valid]
    probs <- .tirt_cat_probs(sp, thj)          # n x n_cat
    scores <- 0:(sp$n_cat - 1)
    E  <- as.vector(probs %*% scores)
    W  <- as.vector(probs %*% (scores^2)) - E^2      # model variance
    W  <- pmax(W, 1e-8)
    # fourth central moment (kurtosis term) C_i = sum_k (k - E)^4 P_k
    cen4 <- outer(E, scores, function(e, s) (s - e)^4)
    C  <- rowSums(probs * cen4)

    resid <- xj - E
    z2 <- resid^2 / W

    outfit <- mean(z2)
    infit  <- sum(resid^2) / sum(W)

    # Wilson-Hilferty standardization
    n <- length(xj)
    q_out2 <- sum(C / W^2) / n^2 - 1 / n
    q_out  <- sqrt(max(q_out2, 1e-8))
    q_in2  <- sum(C - W^2) / (sum(W)^2)
    q_in   <- sqrt(max(q_in2, 1e-8))

    outfit_z <- (outfit^(1/3) - 1) * (3 / q_out) + (q_out / 3)
    infit_z  <- (infit^(1/3)  - 1) * (3 / q_in)  + (q_in  / 3)

    out[[j]] <- data.frame(
      item = parsed$item[j], n = n,
      outfit = round(outfit, 3), outfit_z = round(outfit_z, 3),
      infit = round(infit, 3),  infit_z = round(infit_z, 3),
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, out)
  row.names(res) <- NULL
  res
}
