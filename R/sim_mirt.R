#' Simulate Multidimensional Item Response Theory Data
#'
#' @description
#' Simulates responses from compensatory multidimensional item response models,
#' the multidimensional counterpart of \code{\link{sim_irt}}. It supports the
#' multidimensional binary models estimated by \code{\link{mirt_binary}}
#' (M-Rasch, M-2PL, M-3PL) as well as multidimensional polytomous models
#' (multidimensional GPCM and GRM). Each block of items may load on any subset of
#' the latent dimensions, so both simple-structure and complex-structure designs
#' are easy to generate.
#'
#' @param n_people Integer. Number of examinees.
#' @param item_structure List of lists defining item blocks (see Details).
#' @param dimension Integer. Number of latent dimensions (D).
#' @param theta Numeric matrix (Optional). An \code{n_people} x \code{dimension}
#'   matrix of abilities. If supplied it is used directly.
#' @param theta_mean Numeric. Mean of the latent traits, recycled across
#'   dimensions (used when \code{theta} is \code{NULL}).
#' @param Sigma Numeric matrix (Optional). A \code{dimension} x \code{dimension}
#'   covariance matrix for the latent traits (used when \code{theta} is
#'   \code{NULL}). Defaults to the identity matrix (uncorrelated dimensions).
#'
#' @details
#' Each element of \code{item_structure} is a list describing one block:
#' \itemize{
#'   \item \code{model}: one of \code{"M2PL"}, \code{"M3PL"}, \code{"MRasch"}
#'     (dichotomous) or \code{"MGPCM"}, \code{"MGRM"} (polytomous). The
#'     non-prefixed names (\code{"2PL"}, \code{"GRM"}, ...) are also accepted.
#'   \item \code{n_items}: number of items in the block.
#'   \item \code{dims}: integer vector of the dimensions the block loads on
#'     (default: all dimensions). Slopes on the remaining dimensions are 0.
#'   \item \code{a}: discrimination/slope, given as a range \code{c(lo, hi)} to
#'     sample from, or a single fixed value (default \code{c(0.8, 1.8)}).
#'   \item \code{d}: intercept, given as a range or a fixed value (default
#'     \code{c(-1.5, 1.5)}). Note the compensatory model uses the intercept
#'     metric \eqn{z = a'\theta + d}, matching \code{mirt_binary()}.
#'   \item \code{c}: lower asymptote for \code{"M3PL"} (default \code{0.2}).
#'   \item \code{categories}: number of categories for polytomous blocks
#'     (default \code{3}).
#' }
#'
#' @return A list containing:
#'   \item{resp}{data.frame of responses (rows = people, cols = items).}
#'   \item{true_params}{data.frame of true item parameters, including one
#'     \code{a_Dim} column per dimension, the intercept \code{d}, guessing, and
#'     step/threshold columns for polytomous items.}
#'   \item{theta}{the \code{n_people} x \code{dimension} matrix of true abilities.}
#'
#' @seealso \code{\link{mirt_binary}}, \code{\link{sim_irt}}
#'
#' @examples
#'   # Two correlated dimensions, simple structure
#'   set.seed(2025)
#'   Sigma <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
#'
#'   design <- list(
#'     list(model = "M2PL", n_items = 8, dims = 1),          # loads on Dim 1
#'     list(model = "M2PL", n_items = 8, dims = 2),          # loads on Dim 2
#'     list(model = "M3PL", n_items = 4, dims = c(1, 2), c = 0.2)  # both dims
#'   )
#'
#'   sim <- sim_mirt(n_people = 500, item_structure = design,
#'                   dimension = 2, Sigma = Sigma)
#'
#'   head(sim$resp)
#'   head(sim$true_params)
#'
#'   \donttest{
#'   # Recover with the multidimensional estimator
#'   Q <- as.matrix(sim$true_params[, c("a_Dim1", "a_Dim2")] != 0) * 1
#'   fit <- mirt_binary(sim$resp, model = "2PL", dimension = 2,
#'                      control = list(Q_matrix = Q, max_iter = 10, verbose = FALSE))
#'   head(fit$item_params)
#'   }
#' @export
sim_mirt <- function(n_people = 1000,
                     item_structure = list(),
                     dimension = 2,
                     theta = NULL,
                     theta_mean = 0,
                     Sigma = NULL) {

  # --- 1. Validation ---
  if (!is.numeric(n_people) || n_people < 1) stop("Error: 'n_people' must be a positive integer.")
  if (length(item_structure) == 0) stop("Error: 'item_structure' cannot be empty.")
  D <- as.integer(dimension)
  if (D < 1) stop("Error: 'dimension' must be a positive integer.")

  message("----------------------------------------------------------------")
  message(sprintf("Starting Multidimensional Simulation (N = %d, D = %d)...", n_people, D))

  # --- 2. Generate Abilities ---
  if (!is.null(theta)) {
    theta <- as.matrix(theta)
    if (nrow(theta) != n_people || ncol(theta) != D) {
      stop(sprintf("Error: 'theta' must be a %d x %d matrix.", n_people, D))
    }
    true_theta <- theta
    message(">> Ability (Theta): Using user-supplied matrix.")
  } else {
    if (is.null(Sigma)) Sigma <- diag(D)
    if (!all(dim(Sigma) == c(D, D))) stop(sprintf("Error: 'Sigma' must be a %d x %d matrix.", D, D))
    mu <- rep(theta_mean, length.out = D)
    Z <- matrix(rnorm(n_people * D), n_people, D)
    true_theta <- Z %*% chol(Sigma)
    true_theta <- sweep(true_theta, 2, mu, "+")
    message(sprintf(">> Ability (Theta): Generated from MVN with mean %.2f and supplied Sigma.", theta_mean))
  }

  valid_dich <- c("MRASCH", "M2PL", "M3PL")
  valid_poly <- c("MGPCM", "MGRM")

  # normalize a model label to its multidimensional form
  norm_model <- function(m) {
    mm <- toupper(m)
    switch(mm,
           "RASCH" = , "MRASCH" = "MRASCH",
           "2PL" = , "M2PL" = "M2PL",
           "3PL" = , "M3PL" = "M3PL",
           "GPCM" = , "MGPCM" = "MGPCM",
           "PCM" = "MGPCM",
           "GRM" = , "MGRM" = "MGRM",
           mm)
  }

  response_list <- list()
  param_list <- list()
  current_item_idx <- 1

  # smart parameter parser (mirrors sim_irt)
  get_param <- function(block, name, default, n_req, friendly) {
    val <- block[[name]]
    if (is.null(val)) {
      if (length(default) == 2) return(runif(n_req, default[1], default[2]))
      return(rep(default, n_req))
    }
    if (length(val) == n_req) return(val)
    if (length(val) == 2 && val[1] < val[2]) return(runif(n_req, val[1], val[2]))
    if (length(val) == 1) return(rep(val, n_req))
    stop(sprintf("Parameter '%s' length mismatch (expected %d, or a range).", name, n_req))
  }

  for (i in seq_along(item_structure)) {
    block <- item_structure[[i]]
    if (is.null(block$model)) stop(sprintf("Block %d Error: Missing 'model'.", i))
    if (is.null(block$n_items)) stop(sprintf("Block %d Error: Missing 'n_items'.", i))

    model <- norm_model(block$model)
    n_items <- block$n_items
    if (!model %in% c(valid_dich, valid_poly)) {
      stop(sprintf("Block %d Error: Invalid model '%s'.", i, block$model))
    }

    # dimensions this block loads on
    dims <- block$dims
    if (is.null(dims)) dims <- seq_len(D)
    if (any(dims < 1 | dims > D)) stop(sprintf("Block %d Error: 'dims' out of range.", i))

    # categories
    if (model %in% valid_dich) {
      cats <- 2
    } else {
      cats <- if (is.null(block$categories)) 3 else block$categories
      if (cats < 2) stop(sprintf("Block %d Error: 'categories' must be >= 2.", i))
    }

    message(sprintf(">> Block %d: %d items (%s), loading on dim(s) %s",
                    i, n_items, model, paste(dims, collapse = ",")))

    # slopes: one value per loaded dimension per item
    a_mat <- matrix(0, n_items, D)
    if (model == "MRASCH") {
      a_mat[, dims] <- 1
    } else {
      for (dd in dims) a_mat[, dd] <- get_param(block, "a", c(0.8, 1.8), n_items, "Slope a")
    }

    d_int <- get_param(block, "d", c(-1.5, 1.5), n_items, "Intercept d")
    c_p <- if (model == "M3PL") get_param(block, "c", 0.2, n_items, "Guessing c") else rep(0, n_items)

    block_params <- data.frame(
      item_id = paste0("item_", current_item_idx:(current_item_idx + n_items - 1)),
      block = i,
      model = model,
      categories = cats,
      stringsAsFactors = FALSE
    )
    for (dd in seq_len(D)) block_params[[paste0("a_Dim", dd)]] <- a_mat[, dd]
    block_params$d <- d_int
    block_params$guessing <- c_p

    block_resp <- matrix(NA, n_people, n_items)
    theta_a <- true_theta   # N x D

    if (cats == 2) {
      for (j in seq_len(n_items)) {
        z <- as.vector(theta_a %*% a_mat[j, ]) + d_int[j]
        prob <- c_p[j] + (1 - c_p[j]) / (1 + exp(-z))
        block_resp[, j] <- ifelse(runif(n_people) < prob, 1, 0)
      }
    } else {
      for (j in seq_len(n_items)) {
        comp <- as.vector(theta_a %*% a_mat[j, ])         # composite ability
        d_k <- sort(rnorm(cats - 1, 0, 1), decreasing = TRUE) + d_int[j]
        for (k in seq_len(cats - 1)) block_params[j, paste0("step_", k)] <- d_k[k]

        if (model == "MGRM") {
          prob_cum <- matrix(0, n_people, cats + 1)
          prob_cum[, 1] <- 1
          for (k in seq_len(cats - 1)) prob_cum[, k + 1] <- 1 / (1 + exp(-(comp + d_k[k])))
          probs <- prob_cum[, seq_len(cats)] - prob_cum[, seq_len(cats) + 1]
        } else {                                          # MGPCM
          numer <- matrix(0, n_people, cats)
          running <- 0
          for (k in seq_len(cats - 1)) {
            running <- running + (comp + d_k[k])
            numer[, k + 1] <- running
          }
          numer <- numer - apply(numer, 1, max)
          probs <- exp(numer); probs <- probs / rowSums(probs)
        }
        cum_probs <- t(apply(probs, 1, cumsum))
        rand_vals <- runif(n_people)
        block_resp[, j] <- apply(cum_probs >= rand_vals, 1, function(x) match(TRUE, x)) - 1
      }
    }

    response_list[[i]] <- block_resp
    param_list[[i]] <- block_params
    current_item_idx <- current_item_idx + n_items
  }

  # --- Final Assembly ---
  message("----------------------------------------------------------------")
  message("Constructing final data frames...")

  all_resp <- as.data.frame(do.call(cbind, response_list))
  all_cols <- unique(unlist(lapply(param_list, names)))
  all_params <- do.call(rbind, lapply(param_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) df[missing_cols] <- NA
    df[, all_cols]
  }))
  row.names(all_resp) <- NULL
  row.names(all_params) <- NULL
  colnames(all_resp) <- all_params$item_id
  colnames(true_theta) <- paste0("theta_Dim", seq_len(D))

  message("Simulation Complete.")
  message(sprintf("Summary: %d items, %d examinees, %d dimensions.", ncol(all_resp), n_people, D))
  message("----------------------------------------------------------------")

  list(resp = all_resp, true_params = all_params, theta = true_theta)
}
