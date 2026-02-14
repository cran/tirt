#' Binary (Dichotomous) Item Response Theory Estimation Using Likelihood or Bayesian
#'
#' @description
#' Estimates item and person parameters for binary item response models
#' using either Marginal Maximum Likelihood or Joint Maximum Likelihood.
#' Now supports flexible prior distributions for Bayesian estimation (MAP estimation).
#'
#' @param data A N x J data.frame of dichotomous responses (0/1).
#' @param model String. "Rasch", "2PL" (2-Parameter Logistic), or "3PL" (3-Parameter Logistic).
#' @param method String. "EM" (Marginal Maximum Likelihood via Expectation-Maximization) or "MLE" (Joint Maximum Likelihood). However, using Bayesian will override the likelihood estimation.
#' @param control A \code{list} of control parameters for the estimation algorithm:
#'   \itemize{
#'     \item \code{max_iter}: Maximum number of EM iterations (default = 100).
#'     \item \code{converge_tol}: Convergence criterion for parameter change (default = 1e-4).
#'     \item \code{theta_range}: Numeric vector of length 2 specifying the integration
#'           grid bounds (default = c(-4, 4)).
#'     \item \code{quad_points}: Number of quadrature points (default = 21).
#'     \item \code{verbose}: Logical; if \code{TRUE}, prints progress to console.
#'     \item \code{prior}: A list specifying prior distributions or fixed values for item parameters. Default is NULL (no priors).
#'       For Rasch: \code{list(b = value)} or \code{list(b = function(x) dnorm(x, 0, 1, log=TRUE))}.
#'       For 2PL: \code{list(a = ..., b = ...)}.
#'       For 3PL: \code{list(a = ..., b = ..., g = ...)}.
#'       Each parameter can be:
#'       - A single numeric value (applied to all items)
#'       - A numeric vector of length k (item-specific fixed values)
#'       - A function returning log-density (applied to all items)
#'       - A list of length k with functions or values (item-specific)
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{item_params}: A data frame of estimated item parameters (discrimination, difficulty, guessing) and their standard errors.
#'   \item \code{person_params}: A data frame of estimated person abilities (theta) and standard errors.
#'   \item \code{model_fit}: A data frame containing fit statistics such as Akaike's Information Criterion (AIC), the Bayesian Information Criterion (BIC), and Log-Likelihood.
#'   \item \code{settings}: A list of control parameters used in the estimation.
#' }
#' @examples
#'   # # Simulate data
#'   set.seed(123)
#'   N <- 500; J <- 10
#'   true_theta <- rnorm(N)
#'   true_b <- seq(-2, 2, length.out=J)
#'   true_a <- runif(J, 0.8, 1.2)
#'   data_mat <- matrix(NA, N, J)
#'   for(i in 1:N) {
#'     p <- 1 / (1 + exp(-true_a * (true_theta[i] - true_b)))
#'     data_mat[i,] <- rbinom(J, 1, p)
#'   }
#'   df <- as.data.frame(data_mat)
#'   names(df) <- paste0("Q", 1:J)
#'
#'   # # Run Function without prior
#'   res <- binary_irt(df, model="2PL", method="EM")
#'
#' \donttest{
#'   # # Run Function with prior (function-based)
#'   res_prior <- binary_irt(df, model="2PL", method="EM",
#'                           control=list(prior=list(
#'                             a = function(x) dlnorm(x, 0, 0.5, log=TRUE),
#'                             b = function(x) dnorm(x, 0, 2, log=TRUE)
#'                           )))
#'
#'   # # Run Function with fixed value prior
#'   res_fixed <- binary_irt(df, model="2PL", method="EM",
#'                           control=list(prior=list(
#'                             a = 1,  # Fix all discrimination prior to 1
#'                             b = function(x) dnorm(x, 0, 2, log=TRUE)
#'                           )))
#'
#'   # View Results
#'   head(res$item_params)
#'   head(res$person_params)
#'   print(res$model_fit)
#'   # --- Example 2: With Package Data ---
#'   data("ela1", package = "tirt")
#'   # Subset the first 30 columns (must use the object name 'data_binary')
#'   df <- ela1[, 1:30]
#'   # Run Function on package data
#'   real_res <- binary_irt(df, model="2PL", method="EM")
#'   head(real_res$item_params)
#'   }
#' @export
binary_irt <- function(data,
                       model = "2PL",
                       method = "EM",
                       control = list()) {

  # --- 1. Setup and Defaults ---
  con <- list(
    max_iter = 100,
    converge_tol = 1e-4,
    theta_range = c(-4, 4),
    quad_points = 21,
    nr_max_iter = 20,
    nr_damp = 1.0,
    verbose = TRUE,
    prior = NULL  # No prior by default
  )
  con[names(control)] <- control

  # Extract prior from control
  prior <- con$prior

  # Normalize model name to uppercase
  model <- toupper(model)

  # Define the set of supported models
  valid_models <- c("RASCH", "2PL", "3PL")

  # Validate that every element in the input vector is within the allowed set
  if (!all(model %in% valid_models)) {
    stop("Model must be 'Rasch', '2PL', or '3PL' (case-insensitive)")
  }

  # --- Input Validation ---
  if(!is.data.frame(data)) stop("Input must be a data frame.")

  raw_data <- as.matrix(data)
  if(!is.numeric(raw_data)) {
    stop("Input data must contain only numeric values (0/1). Please remove ID columns or convert non-numeric data.")
  }

  item_names <- colnames(data)
  if(is.null(item_names)) item_names <- paste0("Item_", 1:ncol(data))

  N <- nrow(raw_data)
  J <- ncol(raw_data)

  # --- Prior Validation and Processing ---
  # This function validates and expands priors to match the number of items
  validate_and_expand_prior <- function(prior, model, J, verbose) {
    # If no prior, return NULL
    if(is.null(prior)) {
      if(con$verbose) cat("No prior specified. Using Maximum Likelihood Estimation as Specified.\n")
      return(NULL)
    }

    # Check that prior is a list
    if(!is.list(prior)) {
      stop("Prior must be a list. See documentation for proper format.")
    }

    # 1. Define allowed parameters based on model
    allowed_params <- switch(model,
                             "RASCH" = "b",
                             "2PL"   = c("a", "b"),
                             "3PL"   = c("a", "b", "g"),
                             stop("Unknown model type specified. Support only Rasch, 2PL, and 3PL now..."))

    # 2. Check for extra/unexpected parameters
    extra_params <- setdiff(names(prior), allowed_params)
    if(length(extra_params) > 0) {
      warning(sprintf("Prior specification contains parameter(s) [%s] which are not used in a %s model. These will be ignored.",
                      paste(extra_params, collapse=", "), model))
    }

    # 3. Filter to only allowed and provided parameters
    relevant_params <- intersect(names(prior), allowed_params)
    if(length(relevant_params) == 0) return(NULL)

    if(con$verbose) {
      cat(sprintf("Using MAP estimation for parameters: %s\n", paste(relevant_params, collapse=", ")))
    }

    # Internal check for valid types
    is_valid_prior <- function(x) {
      if(is.numeric(x)) return(TRUE)
      if(is.function(x)) return(TRUE)
      if(is.list(x)) return(all(sapply(x, function(elem) is.numeric(elem) || is.function(elem))))
      return(FALSE)
    }

    expanded_prior <- list()
    for(param in relevant_params) {
      param_prior <- prior[[param]]

      # Validate type (numeric, function, or list)
      isValid <- is.numeric(param_prior) || is.function(param_prior) || is.list(param_prior)
      if(!isValid) stop(sprintf("Prior for '%s' is an invalid type.", param))

      if(!is_valid_prior(param_prior)) {
        stop(sprintf("Prior for parameter '%s' must be numeric, a function, or a list.", param))
      }

      # Case 1: Single numeric value - apply to all items
      if(is.numeric(param_prior) && length(param_prior) == 1) {
        expanded_prior[[param]] <- replicate(J, param_prior, simplify=FALSE)
        if(con$verbose) cat(sprintf("  - Parameter '%s': Fixed value %.3f applied to all %d items\n",
                                    param, param_prior, J))

        # Case 2: Numeric vector of length J - item-specific fixed values
      } else if(is.numeric(param_prior) && length(param_prior) == J) {
        expanded_prior[[param]] <- as.list(param_prior)
        if(con$verbose) cat(sprintf("  - Parameter '%s': Item-specific fixed values for %d items\n",
                                    param, J))

        # Case 3: Numeric vector of wrong length
      } else if(is.numeric(param_prior) && length(param_prior) != J && length(param_prior) != 1) {
        stop(sprintf("Prior for parameter '%s' has length %d, but must be either 1 (applied to all items) or %d (one per item)",
                     param, length(param_prior), J))

        # Case 4: Single function - apply to all items
      } else if(is.function(param_prior)) {
        # Test the function to make sure it works
        test_val <- tryCatch({
          param_prior(0.5)
        }, error = function(e) {
          stop(sprintf("Prior function for parameter '%s' failed when tested with value 0.5. Error: %s",
                       param, e$message))
        })
        if(!is.numeric(test_val) || length(test_val) != 1) {
          stop(sprintf("Prior function for parameter '%s' must return a single numeric value (log-density)", param))
        }

        expanded_prior[[param]] <- replicate(J, param_prior, simplify=FALSE)
        if(con$verbose) cat(sprintf("  - Parameter '%s': Prior function applied to all %d items\n",
                                    param, J))

        # Case 5: List of length J - item-specific priors (can be mix of values and functions)
      } else if(is.list(param_prior) && length(param_prior) == J) {
        # Validate each element
        for(j in 1:J) {
          elem <- param_prior[[j]]
          if(is.function(elem)) {
            # Test the function
            test_val <- tryCatch({
              elem(0.5)
            }, error = function(e) {
              stop(sprintf("Prior function for parameter '%s'[Item %d] failed when tested. Error: %s",
                           param, j, e$message))
            })
            if(!is.numeric(test_val) || length(test_val) != 1) {
              stop(sprintf("Prior function for parameter '%s'[Item %d] must return a single numeric value",
                           param, j))
            }
          } else if(!is.numeric(elem) || length(elem) != 1) {
            stop(sprintf("Prior for parameter '%s'[Item %d] must be a single numeric value or function",
                         param, j))
          }
        }
        expanded_prior[[param]] <- param_prior
        if(con$verbose) cat(sprintf("  - Parameter '%s': Item-specific priors for %d items\n",
                                    param, J))

        # Case 6: List of wrong length
      } else if(is.list(param_prior)) {
        stop(sprintf("Prior list for parameter '%s' has length %d, but must be length %d (one per item)",
                     param, length(param_prior), J))
      }
    }

    return(expanded_prior)
  }

  # Validate and expand the prior
  prior_expanded <- validate_and_expand_prior(con$prior, model, J, con$verbose)

  # Decide on the message based on the RESULT of validation
  if(is.null(prior_expanded)) {
    if(con$verbose) cat("No valid priors for this model. Using Maximum Likelihood Estimation (MLE).\n")
  } else {
    if(con$verbose) cat("Valid priors detected. Using MAP (Maximum A Posteriori) Estimation.\n")
  }

  # --- Helper function to compute log prior density ---
  log_prior_density <- function(value, prior_spec) {
    # If no prior, return 0 (no contribution to log posterior)
    if(is.null(prior_spec)) return(0)

    # If prior is a fixed numeric value, use a very strong normal prior centered at that value
    if(is.numeric(prior_spec)) {
      # Very tight prior: essentially fixes the parameter
      return(dnorm(value, mean=prior_spec, sd=0.001, log=TRUE))
    }

    # If prior is a function, evaluate it
    if(is.function(prior_spec)) {
      result <- tryCatch({
        prior_spec(value)
      }, error = function(e) {
        warning(sprintf("Prior function evaluation failed for value %.3f. Returning -Inf. Error: %s",
                        value, e$message))
        return(-Inf)
      })

      # if the NR algorithm probes slightly outside [0, 1]
      if(is.na(result) || is.nan(result)) return(-Inf)

      if(!is.numeric(result) || length(result) != 1) {
        warning(sprintf("Prior function returned invalid value. Expected single numeric, got %s",
                        class(result)))
        return(-Inf)
      }

      return(result)
    }

    return(0)  # Shouldn't reach here
  }

  # --- Helper function to compute log prior gradient (derivative) using numerical differentiation ---
  log_prior_gradient <- function(value, prior_spec) {
    # If no prior, return 0
    if(is.null(prior_spec)) return(0)

    # If prior is a fixed numeric value, gradient pulls strongly toward that value
    if(is.numeric(prior_spec)) {
      return(-(value - prior_spec) / (0.001^2))
    }

    # If prior is a function, use numerical differentiation
    if(is.function(prior_spec)) {
      h <- 1e-5
      f_plus <- log_prior_density(value + h, prior_spec)
      f_minus <- log_prior_density(value - h, prior_spec)
      return((f_plus - f_minus) / (2 * h))
    }

    return(0)
  }

  # --- Helper function to compute log prior hessian (second derivative) using numerical differentiation ---
  log_prior_hessian <- function(value, prior_spec) {
    # If no prior, return 0
    if(is.null(prior_spec)) return(0)

    # If prior is a fixed numeric value, hessian is constant (very negative = strong pull)
    if(is.numeric(prior_spec)) {
      return(-1 / (0.001^2))
    }

    # If prior is a function, use numerical differentiation
    if(is.function(prior_spec)) {
      h <- 1e-5
      f_0 <- log_prior_density(value, prior_spec)
      f_plus <- log_prior_density(value + h, prior_spec)
      f_minus <- log_prior_density(value - h, prior_spec)
      return((f_plus - 2*f_0 + f_minus) / (h^2))
    }

    return(0)
  }

  # Identify valid items
  item_means <- colMeans(raw_data, na.rm = TRUE)
  valid_items_idx <- which(item_means > 0 & item_means < 1)
  bad_items_idx <- setdiff(1:J, valid_items_idx)

  # Initialize Parameters
  a <- rep(1, J)
  # Calculate proportion correct, avoiding 0 and 1
  p_bounded <- pmin(pmax(item_means, 0.01), 0.99)
  b <- -log(p_bounded / (1 - p_bounded))
  g <- rep(0, J)

  if (model == "3PL") {
    # For 3PL, difficulty must be adjusted for guessing
    # Starting g at 0.15 is safer than 0.2 or 0
    g <- rep(0.15, J)
    # Corrected b calculation for 3PL: b = -log((p - g)/(1 - p)) / a
    b <- -log(pmax(p_bounded - 0.15, 0.01) / (1 - p_bounded))
  } else {
    g <- rep(0, J)
    b <- -log(p_bounded / (1 - p_bounded))
  }

  if (model == "RASCH") a[] <- 1

  # Initialize Theta safely
  raw_scores <- rowMeans(raw_data, na.rm=TRUE)
  raw_scores[is.nan(raw_scores)] <- 0.5
  theta <- as.vector(scale(qlogis(pmin(pmax(raw_scores, 0.01), 0.99))))
  theta[is.nan(theta)] <- 0

  # Quadrature setup for EM
  if (method == "EM") {
    nodes <- seq(-4, 4, length.out = con$quad_points)
    weights <- dnorm(nodes)
    weights <- weights / sum(weights)
  }

  # --- 2. Helper Functions ---

  get_P <- function(th, a_val, b_val, g_val) {
    z <- a_val * (th - b_val)
    z <- pmin(pmax(z, -30), 30)
    p <- g_val + (1 - g_val) * (1 / (1 + exp(-z)))
    return(p)
  }

  # Updated solve_item_nr to include prior information
  solve_item_nr <- function(r_vec, n_vec, theta_vec, cur_a, cur_b, cur_g, mod_type, item_idx) {
    pa <- cur_a; pb <- cur_b; pg <- cur_g

    # Determine which parameters are active for this model
    if(mod_type == "RASCH") active_idx <- 2
    else if(mod_type == "2PL") active_idx <- 1:2
    else active_idx <- 1:3

    # Get priors for this item (if they exist)
    prior_a <- if(!is.null(prior_expanded) && "a" %in% names(prior_expanded)) prior_expanded$a[[item_idx]] else NULL
    prior_b <- if(!is.null(prior_expanded) && "b" %in% names(prior_expanded)) prior_expanded$b[[item_idx]] else NULL
    prior_g <- if(!is.null(prior_expanded) && "g" %in% names(prior_expanded)) prior_expanded$g[[item_idx]] else NULL

    # Log-likelihood function (now includes prior - making it log-posterior)
    ll_func <- function(p_vec) {
      loc_a <- if(1 %in% active_idx) p_vec[1] else 1
      loc_b <- if(mod_type == "RASCH") p_vec[1] else p_vec[2]
      loc_g <- if(3 %in% active_idx) p_vec[3] else 0

      P_tmp <- get_P(theta_vec, loc_a, loc_b, loc_g)
      P_tmp <- pmin(pmax(P_tmp, 1e-9), 1-1e-9)

      # Log-likelihood from data
      log_lik <- sum(r_vec * log(P_tmp) + (n_vec - r_vec) * log(1 - P_tmp))

      # Add log-prior contributions
      log_prior <- 0
      if(!is.null(prior_a) && mod_type != "RASCH") log_prior <- log_prior + log_prior_density(loc_a, prior_a)
      if(!is.null(prior_b))                        log_prior <- log_prior + log_prior_density(loc_b, prior_b)
      if(!is.null(prior_g) && mod_type == "3PL")   log_prior <- log_prior + log_prior_density(loc_g, prior_g)

      return(log_lik + log_prior)
    }

    # Newton-Raphson iterations
    for(iter in 1:con$nr_max_iter) {
      if(mod_type == "RASCH") curr_vec <- c(pb)
      else if(mod_type == "2PL") curr_vec <- c(pa, pb)
      else curr_vec <- c(pa, pb, pg)

      h <- 1e-4
      grad <- numeric(length(curr_vec))
      hess <- matrix(0, length(curr_vec), length(curr_vec))
      f0 <- ll_func(curr_vec)

      # Compute gradient using finite differences
      for(k in 1:length(curr_vec)) {
        tmp <- curr_vec; tmp[k] <- tmp[k] + h
        grad[k] <- (ll_func(tmp) - f0) / h
      }

      # Compute hessian (diagonal only) using finite differences
      for(k in 1:length(curr_vec)) {
        tmp_kk <- curr_vec; tmp_kk[k] <- tmp_kk[k] + h
        tmp_mk <- curr_vec; tmp_mk[k] <- tmp_mk[k] - h
        hess[k,k] <- (ll_func(tmp_kk) - 2*f0 + ll_func(tmp_mk)) / (h^2)
      }

      # Add small regularization for numerical stability
      diag(hess) <- diag(hess) - 1e-5

      try_step <- try(solve(hess, grad), silent = TRUE)
      if(inherits(try_step, "try-error")) break

      step <- try_step * con$nr_damp

      # Check if step is finite before the if statement
      if(any(!is.finite(step))) break

      # Limit the step size to prevent "shooting" to infinity in one jump
      step <- pmin(pmax(step, -2), 2)

      new_vec <- curr_vec - step

      # Update parameters with constraints
      if(mod_type == "RASCH") {
        pb <- max(min(new_vec[1], 15), -15) # Allow extremely low difficulty
      } else if (mod_type == "2PL") {
        pa <- max(min(new_vec[1], 5), 0.01)
        pb <- max(min(new_vec[2], 15), -15)
      } else {
        pa <- max(min(new_vec[1], 5), 0.01)
        pb <- max(min(new_vec[2], 15), -15)
        # Stay strictly within (0.001, 0.5) so priors like Beta don't crash
        pg <- max(min(new_vec[3], 0.5), 0.001)
      }

      # Check convergence
      if(max(abs(step)) < 1e-5) break
    }

    # Compute standard errors
    se_vec <- tryCatch({
      if(mod_type == "RASCH") curr_vec <- c(pb)
      else if(mod_type == "2PL") curr_vec <- c(pa, pb)
      else curr_vec <- c(pa, pb, pg)

      hess_final <- matrix(0, length(curr_vec), length(curr_vec))
      f0 <- ll_func(curr_vec)

      for(k in 1:length(curr_vec)) {
        tmp_kk <- curr_vec; tmp_kk[k] <- tmp_kk[k] + h
        tmp_mk <- curr_vec; tmp_mk[k] <- tmp_mk[k] - h
        hess_final[k,k] <- (ll_func(tmp_kk) - 2*f0 + ll_func(tmp_mk)) / (h^2)
      }

      # Stability: Add a tiny penalty to the diagonal to ensure invertibility
      diag(hess_final) <- diag(hess_final) - 1e-6

      # Check if Hessian is still invalid
      if(any(is.na(hess_final)) || any(diag(hess_final) >= 0)) {
        return(rep(NA, length(curr_vec)))
      }

      sqrt(diag(solve(-hess_final)))
    }, error = function(e) rep(NA, length(curr_vec)), warning = function(w) rep(NA, length(curr_vec)))

    return(list(a=pa, b=pb, g=pg, se=se_vec))
  }

  # --- 3. Main Estimation Loop ---

  if(con$verbose) {
    cat(sprintf("\nStarting %s Estimation using %s algorithm...\n", model, method))
    if(!is.null(prior_expanded)) {
      cat("Bayesian estimation with priors enabled (MAP estimation).\n")
    }
    cat(sprintf("------------------------------------------------\n"))
  }

  is_converged <- FALSE

  for (iter in 1:con$max_iter) {

    # E-STEP
    if (method == "EM") {
      L_terms <- matrix(0, N, con$quad_points)

      for(j in valid_items_idx) {
        if(is.na(a[j]) | is.na(b[j])) next
        P_jq <- get_P(nodes, a[j], b[j], g[j])
        resp <- raw_data[, j]

        term1 <- log(pmax(P_jq, 1e-10))
        term0 <- log(pmax(1 - P_jq, 1e-10))

        r1 <- which(resp == 1)
        r0 <- which(resp == 0)

        if(length(r1)>0) L_terms[r1,] <- L_terms[r1,] + matrix(term1, length(r1), con$quad_points, byrow=TRUE)
        if(length(r0)>0) L_terms[r0,] <- L_terms[r0,] + matrix(term0, length(r0), con$quad_points, byrow=TRUE)
      }

      F_iq <- exp(L_terms) * matrix(weights, N, con$quad_points, byrow=TRUE)
      sum_Fi <- rowSums(F_iq)
      sum_Fi[sum_Fi == 0] <- 1e-10
      Posterior_iq <- F_iq / sum_Fi

    } else {
      # MLE Person Update (no change needed - person parameters don't have priors)
      for(i in 1:N) {
        valid <- !is.na(raw_data[i,]) & (1:J %in% valid_items_idx)
        if(sum(valid)==0) next
        ti <- theta[i]
        for(k in 1:5) {
          p_vec <- get_P(ti, a, b, g)
          diff <- raw_data[i, valid] - p_vec[valid]
          d1 <- sum(a[valid] * diff)
          d2 <- -sum( (a[valid]^2) * p_vec[valid] * (1-p_vec[valid]) )
          if(abs(d2) < 1e-6) break
          step <- d1/d2
          ti <- ti - step
          ti <- max(min(ti, con$theta_range[2]), con$theta_range[1])
          if(abs(step) < 1e-3) break
        }
        theta[i] <- ti
      }
    }

    # M-STEP (now with prior information)
    max_param_change <- 0

    for (j in valid_items_idx) {
      old_params <- c(a[j], b[j], g[j])
      if (method == "EM") {
        r_jq <- colSums(Posterior_iq * ifelse(!is.na(raw_data[,j]) & raw_data[,j]==1, 1, 0))
        n_jq <- colSums(Posterior_iq * ifelse(!is.na(raw_data[,j]), 1, 0))
        # Pass item index to get the correct prior
        res <- solve_item_nr(r_jq, n_jq, nodes, a[j], b[j], g[j], model, j)
      } else {
        valid <- !is.na(raw_data[,j])
        # Pass item index to get the correct prior
        res <- solve_item_nr(raw_data[valid,j], rep(1, sum(valid)), theta[valid], a[j], b[j], g[j], model, j)
      }
      a[j] <- res$a; b[j] <- res$b; g[j] <- res$g
      max_param_change <- max(max_param_change, abs(c(a[j], b[j], g[j]) - old_params), na.rm=TRUE)
    }

    if(con$verbose) cat(paste0("\rIteration ", iter, ": Max Param Change = ", round(max_param_change, 5), "   "))

    # Only allow convergence if parameters moved at least once
    # and the change is small but finite.
    if(iter > 10 && max_param_change < con$converge_tol) {
      if(con$verbose) cat("\n\n>>> Convergence Confirmation: Model Converged!\n")
      is_converged <- TRUE
      break
    }
  }

  # --- User Friendly Reminder: Max Iterations ---
  if(!is_converged && con$verbose) {
    cat(sprintf("\n\n>>> NOTE: Stopped because max Iteration Reached at %d.\n", con$max_iter))
    cat("    Adjust it in control=list(max_iter = ...) to have higher iteration if needed.\n")
  }

  # --- 4. Post-Hoc Estimation & formatting ---

  if(con$verbose) {
    cat("\n------------------------------------------------\n")
    cat("Estimating Final Person Parameters...\n")
    cat(sprintf("Personal Parameter estimated with range [%s, %s].\n", con$theta_range[1], con$theta_range[2]))
    cat("Adjust it in control=list(theta_range=c(..)) if you want a different range.\n")
  }

  if(length(bad_items_idx) > 0) {
    a[bad_items_idx] <- NA; b[bad_items_idx] <- NA; g[bad_items_idx] <- NA
  }

  # Final Person Parameter Estimation
  final_theta <- numeric(N)
  final_theta_se <- numeric(N)
  n_resp <- numeric(N)

  for(i in 1:N) {
    valid <- !is.na(raw_data[i,]) & !is.na(a) & !is.na(b)
    n_resp[i] <- sum(valid)

    if(n_resp[i] == 0) {
      final_theta[i] <- NA; final_theta_se[i] <- NA; next
    }
    score <- sum(raw_data[i, valid])
    if(score == 0) {
      final_theta[i] <- con$theta_range[1]; final_theta_se[i] <- NA; next
    }
    if(score == n_resp[i]) {
      final_theta[i] <- con$theta_range[2]; final_theta_se[i] <- NA; next
    }

    ti <- if(method=="EM") 0 else theta[i]
    for(k in 1:20) {
      p_vec <- get_P(ti, a[valid], b[valid], g[valid])
      u <- raw_data[i, valid]
      num <- (1-g[valid]) * exp(a[valid]*(ti-b[valid]))
      den <- (1 + exp(a[valid]*(ti-b[valid])))^2
      d_P <- a[valid] * num / den
      denom_info <- p_vec * (1-p_vec)
      denom_info[denom_info < 1e-10] <- 1e-10

      score_func <- sum( (u - p_vec) * d_P / denom_info )
      info_func  <- sum( (d_P^2) / denom_info )

      if(is.na(info_func) || info_func < 1e-9) break
      step <- score_func / info_func
      step <- max(min(step, 2), -2)
      ti <- ti + step
      if(abs(step) < 1e-4) break
    }
    final_theta[i] <- ti
    final_theta_se[i] <- if(!is.na(info_func) && info_func > 1e-9) 1/sqrt(info_func) else NA
  }

  if(con$verbose) cat("Person Parameter Estimation Finished.\n")

  # Final Item Statistics with standard errors
  se_a <- rep(NA, J); se_b <- rep(NA, J); se_g <- rep(NA, J)

  for(j in valid_items_idx) {
    if (method == "EM") {
      r_jq <- colSums(Posterior_iq * ifelse(!is.na(raw_data[,j]) & raw_data[,j]==1, 1, 0))
      n_jq <- colSums(Posterior_iq * ifelse(!is.na(raw_data[,j]), 1, 0))
      res <- solve_item_nr(r_jq, n_jq, nodes, a[j], b[j], g[j], model, j)
    } else {
      valid <- !is.na(raw_data[,j])
      res <- solve_item_nr(raw_data[valid,j], rep(1, sum(valid)), theta[valid], a[j], b[j], g[j], model, j)
    }
    if(model == "RASCH") { se_b[j] <- res$se[1] }
    else if (model == "2PL") { se_a[j] <- res$se[1]; se_b[j] <- res$se[2] }
    else { se_a[j] <- res$se[1]; se_b[j] <- res$se[2]; se_g[j] <- res$se[3] }
  }

  # Fit Statistics (including log-prior in log-likelihood for MAP estimation)
  infit <- rep(NA, J); outfit <- rep(NA, J); log_lik <- 0

  for(j in valid_items_idx) {
    valid_p <- !is.na(final_theta) & !is.na(raw_data[,j])
    if(sum(valid_p) == 0) next
    t_use <- final_theta[valid_p]
    x_use <- raw_data[valid_p, j]
    P <- get_P(t_use, a[j], b[j], g[j])
    P_safe <- pmin(pmax(P, 1e-10), 1-1e-10)
    W <- P * (1-P)
    Z <- (x_use - P) / sqrt(W + 1e-10)
    outfit[j] <- sum(Z^2) / length(Z)
    infit[j] <- sum(W * Z^2) / (sum(W) + 1e-10)

    # Log-likelihood from data
    log_lik <- log_lik + sum(x_use * log(P_safe) + (1-x_use) * log(1-P_safe))
  }

  # Add log-prior contribution to log-likelihood if priors are specified
  if(!is.null(prior_expanded)) {
    log_prior_total <- 0
    for(j in valid_items_idx) {
      if(model != "RASCH" && !is.na(a[j])) {
        log_prior_total <- log_prior_total + log_prior_density(a[j], prior_expanded$a[[j]])
      }
      if(!is.na(b[j])) {
        log_prior_total <- log_prior_total + log_prior_density(b[j], prior_expanded$b[[j]])
      }
      if(model == "3PL" && !is.na(g[j])) {
        log_prior_total <- log_prior_total + log_prior_density(g[j], prior_expanded$g[[j]])
      }
    }
    log_lik <- log_lik + log_prior_total
  }

  # Construct Output - SAFELY
  n_count <- as.vector(colSums(!is.na(raw_data)))
  p_val <- as.vector(colMeans(raw_data, na.rm=TRUE))

  if(model == "RASCH") {
    out_items <- data.frame(item=as.vector(item_names),
                            difficulty=as.vector(round(b,3)),
                            difficulty_se=as.vector(round(se_b,3)),
                            number=n_count,
                            pvalue=as.vector(round(p_val,3)),
                            infit=as.vector(round(infit,3)),
                            outfit=as.vector(round(outfit,3)),
                            stringsAsFactors = FALSE)
  } else if (model == "2PL") {
    out_items <- data.frame(item=as.vector(item_names),
                            discrimination=as.vector(round(a,3)),
                            discrimination_se=as.vector(round(se_a,3)),
                            difficulty=as.vector(round(b,3)),
                            difficulty_se=as.vector(round(se_b,3)),
                            number=n_count,
                            pvalue=as.vector(round(p_val,3)),
                            stringsAsFactors = FALSE)
  } else {
    out_items <- data.frame(item=as.vector(item_names),
                            discrimination=as.vector(round(a,3)),
                            discrimination_se=as.vector(round(se_a,3)),
                            difficulty=as.vector(round(b,3)),
                            difficulty_se=as.vector(round(se_b,3)),
                            guess=as.vector(round(g,3)),
                            guess_se=as.vector(round(se_g,3)),
                            number=n_count,
                            pvalue=as.vector(round(p_val,3)),
                            stringsAsFactors = FALSE)
  }

  if(length(bad_items_idx) > 0) {
    cols_to_NA <- which(!names(out_items) %in% c("item", "number", "pvalue"))
    out_items[bad_items_idx, cols_to_NA] <- NA
    if(con$verbose) message("\nWarning: Some items have 0 variance. Parameters set to NA.")
  }

  out_persons <- data.frame(ability=as.vector(round(final_theta,3)),
                            ability_se=as.vector(round(final_theta_se,3)),
                            number=as.vector(n_resp),
                            stringsAsFactors = FALSE)

  k_params <- length(valid_items_idx) * (if(model=="RASCH") 1 else if(model=="2PL") 2 else 3)
  aic <- 2*k_params - 2*log_lik
  bic <- k_params*log(N) - 2*log_lik

  # Add note about MAP estimation in fit statistics
  log_lik_label <- if(!is.null(prior_expanded)) "LogPosterior" else "LogLikelihood"

  out_fit <- data.frame(Index=c(log_lik_label, "AIC", "BIC", "Iterations"),
                        Value=as.vector(c(round(log_lik,3), round(aic,3), round(bic,3), iter)),
                        stringsAsFactors = FALSE)

  # Settings summary for user reference
  estimation_type <- if(!is.null(prior_expanded)) "MAP (with priors)" else "MLE (no priors)"
  out_settings <- data.frame(
    Parameter = c("Model", "Method", "Estimation", "Max_Iter", "Converge_Tol", "Theta_Min", "Theta_Max"),
    Value = as.vector(c(model, method, estimation_type, con$max_iter, con$converge_tol,
                        con$theta_range[1], con$theta_range[2])),
    stringsAsFactors = FALSE
  )

  if(con$verbose) {
    cat("Finished All Estimation.\n")
    cat("------------------------------------------------\n")
  }

  return(list(item_params=out_items,
              person_params=out_persons,
              model_fit=out_fit,
              settings=out_settings))
}
