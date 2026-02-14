#' Polytomous Item Response Theory Estimation Using Likelihood or Bayesian
#'
#' @description
#' Estimates item and person parameters for polytomous item response theory models
#' using either Marginal Maximum Likelihood or Joint Maximum Likelihood.
#' Now supports flexible prior distributions for Bayesian estimation (MAP estimation).
#'
#' @param data A N x J data.frame of polytomous responses (0, 1, 2...).
#'             Missing values should be NA. Categories must be continuous integers.
#' @param model String. "GPCM" (Generalized Partial Credit Model), "PCM" (Partial Credit Model), or "GRM" (Graded Response Model).
#' @param method String. "EM" (Marginal Maximum Likelihood via Expectation-Maximization) or "MLE" (Joint Maximum Likelihood). However, using Bayesian will override the likelihood estimation.
#' @param control A \code{list} of control parameters for the estimation algorithm:
#'   \itemize{
#'     \item \code{max_iter}: Maximum number of EM iterations (default = 100).
#'     \item \code{converge_tol}: Convergence criterion for parameter change (default = 1e-4).
#'     \item \code{theta_range}: Numeric vector of length 2 specifying the integration
#'           grid bounds (default = c(-4, 4)).
#'     \item \code{quad_points}: Number of quadrature points (default = 21).
#'     \item \code{verbose}: Logical; if \code{TRUE}, prints progress to console.
#'     \item \code{prior}: A list specifying prior distributions for item parameters. Default is NULL (no priors).
#'       For GRM: \code{list(a = function(x) dlnorm(x, 0, 0.5, log=TRUE), d = function(x) dnorm(x, 0, 2, log=TRUE))}.
#'       For GPCM: \code{list(a = function(x) dlnorm(x, 0, 0.5, log=TRUE), d = function(x) dnorm(x, 0, 2, log=TRUE))}.
#'       For PCM: \code{list(d = function(x) dnorm(x, 0, 2, log=TRUE))}.
#'       Each prior is a function returning log-density and applies to ALL items.
#'       Fixed value priors are NOT supported.
#'       Item-specific priors are NOT supported.
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{item_params}: Data frame of estimated parameters (a, thresholds).
#'   \item \code{person_params}: A data frame of estimated person abilities (theta) and standard errors.
#'   \item \code{model_fit}: A data frame containing fit statistics such as Akaike’s Information Criterion (AIC) and the Bayesian Information Criterion (BIC).
#'   \item \code{settings}: A list of control parameters used in the estimation.
#' }
#' @examples
#'   # --- Example 1: Simulation (GPCM) ---
#'   set.seed(2026)
#'   N <- 500; J <- 5
#'   n_cats <- c(3, 4, 3, 5, 4)
#'
#'   true_theta <- rnorm(N)
#'   true_a <- runif(J, 0.8, 1.2)
#'   true_d <- list()
#'
#'   # Generate Thresholds
#'   for(j in 1:J) {
#'     steps <- sort(rnorm(n_cats[j]-1, mean = 0, sd = 1.0))
#'     true_d[[j]] <- c(0, cumsum(steps))
#'   }
#'
#'   # Simulation Helper (GPCM Logic)
#'   generate_resp <- function(theta, a, d_vec, n_cat) {
#'     probs <- matrix(0, length(theta), n_cat)
#'     for(k in 1:n_cat) {
#'       z <- a * (k-1) * theta - d_vec[k]
#'       probs[,k] <- exp(z)
#'     }
#'     probs <- probs / rowSums(probs)
#'     apply(probs, 1, function(p) sample(0:(n_cat-1), 1, prob=p))
#'   }
#'
#'   # Create Data
#'   sim_data <- matrix(NA, nrow = N, ncol = J)
#'   for(j in 1:J) {
#'     sim_data[,j] <- generate_resp(true_theta, true_a[j], true_d[[j]], n_cats[j])
#'   }
#'   df_sim <- as.data.frame(sim_data)
#'
#'   # Run Estimation (GPCM to match simulation logic without prior)
#'   res <- polytomous_irt(df_sim, model="GPCM", method="EM",
#'                         control=list(max_iter=20, verbose=TRUE))
#'
#'   head(res$item_params)
#'   print(res$model_fit)
#'
#'   \donttest{
#'   # Run Estimation with prior (MAP)
#'   res_prior <- polytomous_irt(df_sim, model="PCM", method="EM",
#'                               control=list(max_iter=20, verbose=FALSE,
#'                                           prior=list(
#'                                             d = function(x) dnorm(x, 0, 2, log=TRUE)
#'                                           )))
#'   head(res$item_params)
#'   print(res$model_fit)
#'
#'   # --- Example 2: With Package Data (GRM) ---
#'   data("ela1", package = "tirt")
#'
#'   # Subset polytomous items (columns 31 to 45)
#'   df_poly <- ela1[, 31:45]
#'
#'   # Run Estimation using GRM
#'   real_res <- polytomous_irt(df_poly, model="GRM", method="EM",
#'                              control = list(max_iter = 1000))
#'
#'   head(real_res$item_params)
#'   head(real_res$person_params)
#'   print(real_res$model_fit)
#'   # Run Estimation using GRM with prior
#'   real_res2 <- polytomous_irt(df_poly, model="GRM", method="EM",
#'                              control = list(max_iter = 1000,
#'                                            prior = list(
#'                                              a = function(x) dlnorm(x, 0, 0.5, log=TRUE),
#'                                              d = function(x) dnorm(x, 0, 2, log=TRUE)
#'                                            )))
#'
#'   head(real_res2$item_params)
#'   head(real_res2$person_params)
#'   print(real_res2$model_fit)
#'   }
#' @export
polytomous_irt <- function(data,
                           model = "GPCM",
                           method = "EM",
                           control = list()) {

  # --- 1. Setup and Defaults ---
  con <- list(
    max_iter = 100,
    converge_tol = 1e-4,
    theta_range = c(-4, 4),
    quad_points = 21,
    nr_max_iter = 20,
    nr_damp = 0.8,
    verbose = TRUE
  )
  con[names(control)] <- control

  # Extract prior from control
  prior <- con$prior

  # Normalize model name to uppercase
  model <- toupper(model)

  # Define the set of supported models
  valid_models <- c("GPCM", "PCM", "GRM")

  # Validate that every element in the input vector is within the allowed set
  if (!all(model %in% valid_models)) {
    stop("Model must be 'GPCM', 'PCM', or 'GRM' (case-insensitive)")
  }

  # --- Input Validation ---
  if(!is.data.frame(data)) stop("Input must be a data frame.")

  # Check data types
  raw_mat <- as.matrix(data)
  if(!is.numeric(raw_mat)) stop("Data must be numeric. Remove ID columns.")

  item_names <- colnames(data)
  if(is.null(item_names)) item_names <- paste0("Item_", 1:ncol(data))

  N <- nrow(raw_mat)
  J <- ncol(raw_mat)

  # --- Data Preprocessing (Shift to 0-based) ---
  min_val <- min(raw_mat, na.rm=TRUE)
  if(is.finite(min_val) && min_val > 0) {
    if(con$verbose) cat("Note: Responses shifted to start at 0.\n")
    raw_mat <- raw_mat - min_val
  }

  # Determine max categories per item
  n_cats <- apply(raw_mat, 2, function(x) {
    if(all(is.na(x))) return(0)
    max(x, na.rm=TRUE) + 1
  })
  max_K <- max(n_cats)
  if(max_K < 2) stop("Data must contain at least 2 categories.")

  # Identify valid items (Variance > 0)
  item_vars <- apply(raw_mat, 2, var, na.rm=TRUE)
  valid_items_idx <- which(item_vars > 0)

  if(length(valid_items_idx) == 0) stop("No valid items found.")

  # --- Prior Validation and Processing ---
  # This function validates priors for polytomous models
  validate_and_expand_prior <- function(prior, model, verbose) {
    # If no prior, return NULL
    if(is.null(prior)) {
      if(verbose) cat("No prior specified. Using Maximum Likelihood Estimation as Specified.\n")
      return(NULL)
    }

    # Check that prior is a list
    if(!is.list(prior)) {
      stop("Prior must be a list. See documentation for proper format.")
    }

    # 1. Define allowed parameters based on model
    allowed_params <- switch(model,
                             "PCM"  = "d",
                             "GPCM" = c("a", "d"),
                             "GRM"  = c("a", "d"),
                             stop("Unknown model type specified. Support only PCM, GPCM, and GRM now..."))

    # 2. Check for extra/unexpected parameters
    extra_params <- setdiff(names(prior), allowed_params)
    if(length(extra_params) > 0) {
      warning(sprintf("Prior specification contains parameter(s) [%s] which are not used in a %s model. These will be ignored.",
                      paste(extra_params, collapse=", "), model))
    }

    # 3. Filter to only allowed and provided parameters
    relevant_params <- intersect(names(prior), allowed_params)
    if(length(relevant_params) == 0) return(NULL)

    if(verbose) {
      cat(sprintf("Using MAP (Bayesian) estimation for parameters: %s\n", paste(relevant_params, collapse=", ")))
    }

    # 4. Validate each prior is a function (no fixed values, no item-specific)
    validated_prior <- list()
    for(param in relevant_params) {
      param_prior <- prior[[param]]

      # Check if it's a function
      if(!is.function(param_prior)) {
        stop(sprintf("Prior for parameter '%s' must be a function returning log-density. Fixed value priors and item-specific priors are not supported for polytomous models.",
                     param))
      }

      # Test the function
      test_val <- tryCatch({
        param_prior(0.5)
      }, error = function(e) {
        stop(sprintf("Prior function for parameter '%s' failed when tested: %s", param, e$message))
      })

      if(!is.numeric(test_val) || length(test_val) != 1) {
        stop(sprintf("Prior function for parameter '%s' must return a single numeric value (log-density).", param))
      }

      validated_prior[[param]] <- param_prior

      if(verbose) {
        cat(sprintf("  - Parameter '%s': Prior function applied to ALL items/parameters\n", param))
      }
    }

    return(validated_prior)
  }
  # Validate and process priors
  prior_expanded <- validate_and_expand_prior(prior, model, con$verbose)

  # Helper function to compute log-prior density
  log_prior_density <- function(value, prior_func) {
    if(is.null(prior_func)) return(0)
    if(!is.function(prior_func)) return(0)
    if(is.na(value)) return(0)
    tryCatch(prior_func(value), error = function(e) 0)
  }

  # --- Initialize Parameters ---
  a <- rep(1, J)
  d <- matrix(0, J, max_K)

  # --- UPDATED: Initialize Standard Error Storage ---
  se_a <- rep(NA, J)
  se_d <- matrix(NA, J, max_K)

  for(j in valid_items_idx) {
    valid_resp <- raw_mat[!is.na(raw_mat[,j]), j]
    if(length(valid_resp) == 0) next

    props <- table(factor(valid_resp, levels=0:(n_cats[j]-1))) / length(valid_resp)
    cum_props <- cumsum(props)
    cum_props <- pmin(pmax(cum_props, 0.01), 0.99)
    thresh <- -qlogis(cum_props[-length(cum_props)])

    if(model == "GRM") {
      d[j, 1:(n_cats[j]-1)] <- sort(thresh)
    } else {
      K <- n_cats[j]
      steps <- seq(-0.5, 0.5, length.out = K-1)
      d[j, 2:K] <- steps
    }
  }

  if (model == "PCM") a[] <- 1

  # Initialize Theta
  raw_scores <- rowSums(raw_mat, na.rm=TRUE)
  theta <- as.vector(scale(raw_scores))
  theta[is.nan(theta) | is.na(theta)] <- 0

  # Quadrature setup (Only used for EM)
  nodes <- seq(con$theta_range[1], con$theta_range[2], length.out = con$quad_points)
  weights <- dnorm(nodes)
  weights <- weights / sum(weights)


  # --- 2. Helper Functions ---

  get_P_poly <- function(th_vec, a_val, d_vec, mod_type, n_cat) {
    n_q <- length(th_vec)
    probs <- matrix(0, n_q, n_cat)

    if (mod_type %in% c("GPCM", "PCM")) {
      for(k in 1:n_cat) {
        z <- a_val * (k-1) * th_vec - d_vec[k]
        z <- pmin(pmax(z, -50), 50)
        probs[,k] <- exp(z)
      }
      row_sums <- rowSums(probs)
      row_sums[row_sums == 0] <- 1
      probs <- probs / row_sums

    } else if (mod_type == "GRM") {
      P_star <- matrix(0, n_q, n_cat + 1)
      P_star[,1] <- 1
      P_star[,n_cat + 1] <- 0

      eff_d <- d_vec[1:(n_cat-1)]
      if(is.unsorted(eff_d)) eff_d <- sort(eff_d)

      for(k in 1:(n_cat-1)) {
        z <- a_val * (th_vec - eff_d[k])
        z <- pmin(pmax(z, -30), 30)
        P_star[, k+1] <- 1 / (1 + exp(-z))
      }
      for(k in 1:n_cat) {
        probs[,k] <- P_star[,k] - P_star[,k+1]
      }
      probs <- pmax(probs, 1e-10)
    }
    return(probs)
  }

  solve_item_poly <- function(r_mat, n_vec, th_vec, cur_a, cur_d, mod_type, n_cat) {
    if(mod_type == "PCM") {
      param_vec <- cur_d[2:n_cat]
    } else if(mod_type == "GRM") {
      param_vec <- c(cur_a, cur_d[1:(n_cat-1)])
    } else {
      param_vec <- c(cur_a, cur_d[2:n_cat])
    }

    n_params <- length(param_vec)

    ll_func <- function(p_v) {
      if(mod_type == "PCM") {
        l_a <- 1; l_d <- numeric(max_K); l_d[2:n_cat] <- p_v
      } else {
        l_a <- p_v[1]; l_d <- numeric(max_K)
        if(mod_type == "GRM") l_d[1:(n_cat-1)] <- p_v[-1]
        else l_d[2:n_cat] <- p_v[-1]
      }
      probs <- get_P_poly(th_vec, l_a, l_d, mod_type, n_cat)
      probs <- pmin(pmax(probs, 1e-10), 1-1e-10)

      log_lik_data <- sum(r_mat * log(probs))

      # Add priors if specified
      if(!is.null(prior_expanded)) {
        if(mod_type != "PCM" && !is.null(prior_expanded$a)) {
          log_lik_data <- log_lik_data + log_prior_density(l_a, prior_expanded$a)
        }
        if(!is.null(prior_expanded$d)) {
          if(mod_type == "GRM") {
            for(k in 1:(n_cat-1)) {
              log_lik_data <- log_lik_data + log_prior_density(l_d[k], prior_expanded$d)
            }
          } else {
            for(k in 2:n_cat) {
              log_lik_data <- log_lik_data + log_prior_density(l_d[k], prior_expanded$d)
            }
          }
        }
      }

      return(log_lik_data)
    }

    for(iter in 1:con$nr_max_iter) {
      f0 <- ll_func(param_vec)
      h <- 1e-4
      grad <- numeric(n_params)
      hess <- matrix(0, n_params, n_params)

      for(k in 1:n_params) {
        p_up <- param_vec; p_up[k] <- p_up[k] + h
        grad[k] <- (ll_func(p_up) - f0) / h
      }

      for(k in 1:n_params) {
        p_up <- param_vec; p_up[k] <- p_up[k] + h
        p_dn <- param_vec; p_dn[k] <- p_dn[k] - h
        hess[k,k] <- (ll_func(p_up) - 2*f0 + ll_func(p_dn)) / (h^2)
      }
      diag(hess) <- diag(hess) - 1e-4

      try_step <- try(solve(hess, grad), silent = TRUE)
      if(inherits(try_step, "try-error")) break

      step <- try_step * con$nr_damp
      step <- pmax(pmin(step, 1.0), -1.0)
      param_vec <- param_vec - step

      if(mod_type != "PCM") param_vec[1] <- max(param_vec[1], 0.05)
      if(max(abs(step)) < 1e-4) break
    }

    se_vec <- tryCatch({
      suppressWarnings(sqrt(diag(solve(-hess))))
    }, error = function(e) rep(NA, n_params))

    ret_a <- if(mod_type=="PCM") 1 else param_vec[1]
    ret_d <- cur_d
    ret_se_a <- NA
    ret_se_d <- rep(NA, max_K)

    if(mod_type == "PCM") {
      ret_d[2:n_cat] <- param_vec; ret_se_d[2:n_cat] <- se_vec
    } else if (mod_type == "GRM") {
      raw_thresh <- param_vec[-1]
      ret_d[1:(n_cat-1)] <- sort(raw_thresh)
      ret_se_a <- se_vec[1]
      if(is.unsorted(raw_thresh)) {
        ret_se_d[1:(n_cat-1)] <- NA
      } else {
        ret_se_d[1:(n_cat-1)] <- se_vec[-1]
      }
    } else {
      ret_d[2:n_cat] <- param_vec[-1]; ret_se_a <- se_vec[1]; ret_se_d[2:n_cat] <- se_vec[-1]
    }
    return(list(a=ret_a, d=ret_d, se_a=ret_se_a, se_d=ret_se_d))
  }

  # --- 3. Main Estimation Loop ---

  if(con$verbose) cat(sprintf("\nStarting Polytomous Estimation (%s) using %s...\n", model, method))

  converged <- FALSE

  for (iter in 1:con$max_iter) {

    # --- E-Step (EM) or Person Step (MLE) ---
    if (method == "EM") {
      # Calculate Posterior probabilities for Quadrature nodes
      L_mat <- matrix(0, N, con$quad_points)
      for(j in valid_items_idx) {
        P_jq <- get_P_poly(nodes, a[j], d[j,], model, n_cats[j])
        resp <- raw_mat[, j]
        valid_rows <- !is.na(resp)
        p_indices <- resp[valid_rows] + 1
        item_loglik <- matrix(0, sum(valid_rows), con$quad_points)
        for(k in 1:n_cats[j]) {
          is_cat <- which(p_indices == k)
          if(length(is_cat) > 0) {
            val <- log(pmax(P_jq[,k], 1e-50))
            item_loglik[is_cat, ] <- matrix(val, length(is_cat), con$quad_points, byrow=TRUE)
          }
        }
        L_mat[valid_rows, ] <- L_mat[valid_rows, ] + item_loglik
      }
      F_iq <- exp(L_mat) * matrix(weights, N, con$quad_points, byrow=TRUE)
      sum_Fi <- rowSums(F_iq)
      sum_Fi[sum_Fi == 0] <- 1e-10
      Posterior_iq <- F_iq / sum_Fi

    } else {
      # Method == "MLE": Estimate Theta for every person given current items
      for(i in 1:N) {
        valid <- !is.na(raw_mat[i,]) & (1:J %in% valid_items_idx)
        if(sum(valid)==0) next
        ti <- theta[i]
        for(k in 1:5) { # Few NR steps for efficiency
          num <- 0; den <- 0
          for(j in which(valid)) {
            probs <- get_P_poly(ti, a[j], d[j,], model, n_cats[j])[1,]
            x <- raw_mat[i,j]
            h <- 1e-4
            # Numerical derivatives
            p0 <- probs[x+1]
            pp <- get_P_poly(ti+h, a[j], d[j,], model, n_cats[j])[1, x+1]
            pm <- get_P_poly(ti-h, a[j], d[j,], model, n_cats[j])[1, x+1]
            d1 <- (log(pp)-log(pm))/(2*h)
            d2 <- (log(pp)-2*log(p0)+log(pm))/(h^2)
            num <- num + d1; den <- den + d2
          }
          if(is.na(den) || abs(den) < 1e-5) break
          step <- num/den
          ti <- ti - step
          ti <- max(min(ti, con$theta_range[2]), con$theta_range[1])
          if(abs(step) < 1e-3) break
        }
        theta[i] <- ti
      }
      # Re-center Theta to fix scale indeterminacy in JML
      theta <- theta - mean(theta, na.rm=TRUE)
    }

    # --- M-Step (EM) or Item Step (MLE) ---
    max_change <- 0

    for(j in valid_items_idx) {
      old_par <- c(a[j], d[j,])

      if(method == "EM") {
        # --- EM Preparation ---
        # R_qk: Expected number of people in category k at quadrature node q
        R_qk <- matrix(0, con$quad_points, n_cats[j])
        resp <- raw_mat[, j]
        valid <- !is.na(resp)
        for(k in 1:n_cats[j]) {
          idx <- which(valid & resp == (k-1))
          if(length(idx) > 0) R_qk[,k] <- colSums(Posterior_iq[idx, , drop=FALSE])
        }
        N_q <- rowSums(R_qk)
        use_nodes <- nodes # Use global quadrature nodes

      } else {
        # --- MLE (JML) Preparation ---
        valid_idx <- which(!is.na(raw_mat[, j]))
        if(length(valid_idx) == 0) next

        y <- raw_mat[valid_idx, j]
        use_nodes <- theta[valid_idx]

        n_resp <- length(y)
        R_qk <- matrix(0, n_resp, n_cats[j])

        idx_mat <- cbind(1:n_resp, y + 1)
        R_qk[idx_mat] <- 1

        N_q <- rep(1, n_resp)
      }

      # Solve Item Parameters
      res <- solve_item_poly(R_qk, N_q, use_nodes, a[j], d[j,], model, n_cats[j])
      a[j] <- res$a
      d[j,] <- res$d

      # --- UPDATED: Store Standard Errors ---
      se_a[j] <- res$se_a
      se_d[j,] <- res$se_d

      max_change <- max(max_change, abs(c(a[j], d[j,]) - old_par), na.rm=TRUE)
    }

    if(con$verbose) cat(paste0("\rIteration ", iter, ": Max Param Change = ", round(max_change, 5), "   "))

    if(max_change < con$converge_tol) {
      if(con$verbose) cat("\n\n>>> Convergence Confirmation: Model Converged!\n")
      converged <- TRUE
      break
    }
  }

  if(!converged && con$verbose) {
    cat(sprintf("\nNOTE: Stopped because max Iteration Reached at %d.\n", con$max_iter))
    cat("    Adjust it in control=list(max_iter = ...) to have higher iteration if needed.\n")
  }

  # --- 4. Post-Hoc Estimation & Output ---
  if(con$verbose) {
    cat("\n------------------------------------------------\n")
    cat("Estimating Final Person Parameters...\n")
  }

  final_theta <- numeric(N)
  final_se <- numeric(N)

  for(i in 1:N) {
    valid <- !is.na(raw_mat[i,]) & (1:J %in% valid_items_idx)
    if(sum(valid)==0) { final_theta[i] <- NA; final_se[i] <- NA; next }

    # Use existing estimates as starting values for speed
    ti <- if(method=="EM") 0 else theta[i]

    for(k in 1:20) {
      score_num <- 0; info_den <- 0
      for(j in which(valid)) {
        probs <- get_P_poly(ti, a[j], d[j,], model, n_cats[j])[1,]
        x <- raw_mat[i, j]
        h <- 1e-4
        probs_p <- get_P_poly(ti+h, a[j], d[j,], model, n_cats[j])[1,]
        probs_m <- get_P_poly(ti-h, a[j], d[j,], model, n_cats[j])[1,]

        lp0 <- log(max(probs[x+1], 1e-50))
        lpp <- log(max(probs_p[x+1], 1e-50))
        lpm <- log(max(probs_m[x+1], 1e-50))

        d1 <- (lpp - lpm) / (2*h)
        d2 <- (lpp - 2*lp0 + lpm) / (h^2)
        score_num <- score_num + d1; info_den <- info_den - d2
      }
      if(is.na(info_den) || abs(info_den) < 1e-6) break
      step <- score_num / info_den
      step <- max(min(step, 1.0), -1.0)
      ti <- ti + step
      if(abs(step) < 1e-4) break
    }
    final_theta[i] <- max(min(ti, con$theta_range[2]), con$theta_range[1])
    if(!is.na(info_den) && info_den > 1e-6) final_se[i] <- 1/sqrt(info_den) else final_se[i] <- NA
  }

  if(con$verbose) cat("Person Parameter Estimation Finished.\n")

  # --- Output Construction (UPDATED) ---
  out_list <- list(
    item = as.character(item_names),
    discrimination = as.numeric(round(a, 3)),
    discrimination_se = as.numeric(round(se_a, 3)), # <--- _se suffix
    n_cat = as.integer(n_cats)
  )

  n_thresh_cols <- max_K - 1
  if(n_thresh_cols > 0) {
    for(k in 1:n_thresh_cols) {
      col_vals <- rep(NA_real_, J)
      se_col_vals <- rep(NA_real_, J)

      for(j in valid_items_idx) {
        if(model == "GRM") {
          if(k < n_cats[j]) {
            col_vals[j] <- d[j, k]
            se_col_vals[j] <- se_d[j, k]
          }
        } else {
          if(k < n_cats[j]) {
            col_vals[j] <- d[j, k+1]
            se_col_vals[j] <- se_d[j, k+1]
          }
        }
      }
      col_name <- paste0(if(model=="GRM") "thresh_" else "step_", k)
      out_list[[col_name]] <- round(col_vals, 3)

      se_col_name <- paste0(col_name, "_se") # <--- _se suffix
      out_list[[se_col_name]] <- round(se_col_vals, 3)
    }
  }

  # Ensure no rownames in output data.frames
  df_items <- data.frame(out_list, stringsAsFactors = FALSE, row.names = NULL)

  df_persons <- data.frame(
    ability = as.numeric(round(final_theta, 3)),
    ability_se = as.numeric(round(final_se, 3)), # <--- _se suffix
    n_resp = as.integer(rowSums(!is.na(raw_mat))),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  log_lik <- 0
  for(i in 1:N) {
    valid <- !is.na(raw_mat[i,]) & !is.na(final_theta[i])
    if(!any(valid)) next
    for(j in which(valid)) {
      probs <- get_P_poly(final_theta[i], a[j], d[j,], model, n_cats[j])[1,]
      log_lik <- log_lik + log(max(probs[raw_mat[i,j]+1], 1e-50))
    }
  }

  n_par <- length(valid_items_idx) + sum(n_cats[valid_items_idx]-1)
  aic <- 2*n_par - 2*log_lik
  bic <- n_par*log(N) - 2*log_lik

  df_fit <- data.frame(
    Index = c("LogLikelihood", "AIC", "BIC"),
    Value = as.numeric(c(round(log_lik,3), round(aic,3), round(bic,3))),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if(con$verbose) {
    cat("Finished All Estimation.\n")
    cat("------------------------------------------------\n")
  }

  return(list(
    item_params = df_items,
    person_params = df_persons,
    model_fit = df_fit,
    settings = list(model=model, method=method)
  ))
}
