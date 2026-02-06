#' Unidimensional Binary (Dichotomous) Testlet Response Theory Estimation
#'
#' @description
#' Estimates item and person parameters for Unidimensional Binary (Dichotomous) Testlet response models
#' using Penalized Expectation-Maximization or
#' Joint Maximum Likelihood Estimation with stabilization.
#'
#' @param data A data.frame of binary responses (0/1). Rows=persons, Cols=items in testlets.
#' @param group A list defining testlet structures. Example: `list(c(1,2,3), c(4,5,6))`.
#' @param model Character. One of "RaschT" (Rasch Testlet), "2PLT" (2-Parameter Logistic Testlet), "3PLT" (3-Parameter Logistic Testlet), "BiFT" (Bifactor).
#' @param method Character. "EM" (Marginal Maximum Likelihood via Expectation-Maximization) or "MLE" (Joint Maximum Likelihood).
#' @param control A \code{list} of control parameters for the estimation algorithm:
#'   \itemize{
#'     \item \code{max_iter}: Maximum number of EM iterations (default = 100).
#'     \item \code{converge_tol}: Convergence criterion for parameter change (default = 1e-4).
#'     \item \code{theta_range}: Numeric vector of length 2 specifying the integration
#'           grid bounds (default = c(-4, 4)).
#'     \item \code{quad_points}: Number of quadrature points (default = 21).
#'     \item \code{verbose}: Logical; if \code{TRUE}, prints progress to console.
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{item_params}: Estimated item parameters.
#'   \item \code{person_params}: Estimated person abilities and testlet effects.
#'   \item \code{model_fit}: A data frame containing iterations and fit statistics such as Akaike’s Information Criterion (AIC), the Bayesian Information Criterion (BIC), and Log-Likelihood.
#' }
#'
#' @examples
#'   # --- Example: Simulation (2PLT) ---
#'   set.seed(2025)
#'   n_persons <- 500
#'   n_testlets <- 3
#'   items_per_testlet <- 3
#'   n_items <- n_testlets * items_per_testlet
#'
#'   # 1. Generate Parameters
#'   # Discrimination (a): Varying -> 2PLT
#'   a_true <- runif(n_items, 0.8, 1.5)
#'   # Difficulty (b)
#'   b_true <- seq(-1, 1, length.out = n_items)
#'   # Testlet Variances (Sigma)
#'   sigma_true <- c(1.0, 1.5, 2.0)
#'
#'   # 2. Generate Person Params
#'   theta_true <- rnorm(n_persons, 0, 1)
#'   gamma_matrix <- matrix(0, nrow = n_persons, ncol = n_testlets)
#'   for(d in 1:n_testlets) {
#'     gamma_matrix[, d] <- rnorm(n_persons, 0, sigma_true[d])
#'   }
#'
#'   # 3. Generate Responses
#'   resp_matrix <- matrix(0, nrow = n_persons, ncol = n_items)
#'   colnames(resp_matrix) <- paste0("Item_", 1:n_items)
#'   group_list <- list()
#'
#'   idx_counter <- 1
#'   for(d in 1:n_testlets) {
#'     indices <- idx_counter:(idx_counter + items_per_testlet - 1)
#'     group_list[[d]] <- indices
#'
#'     for(i in indices) {
#'       # 2PLT Model: a * (theta + gamma - b)
#'       lin <- a_true[i] * (theta_true + gamma_matrix[, d] - b_true[i])
#'       prob <- 1 / (1 + exp(-lin))
#'       resp_matrix[, i] <- rbinom(n_persons, 1, prob)
#'     }
#'     idx_counter <- idx_counter + items_per_testlet
#'   }
#'   df_sim <- as.data.frame(resp_matrix)
#'
#'   # 4. Run Estimation
#'   # We use "2PLT" because data was generated with varying 'a'
#'   res <- trt_binary(
#'     data = df_sim,
#'     group = group_list,
#'     model = "2PLT",
#'     method = "EM",
#'     control = list(max_iter = 20, verbose = FALSE)
#'   )
#'
#'   head(res$item_params)
#'   head(res$person_params)
#' @export
trt_binary <- function(
    data,
    group,
    model = c("RaschT", "2PLT", "3PLT", "BiFT"),
    method = c("EM", "MLE"),
    control = list()
) {

  # ===========================================================================
  # 1. SETUP & VALIDATION
  # ===========================================================================
  model  <- match.arg(model)
  method <- match.arg(method)

  # Default controls with added regularization parameters
  con <- list(
    max_iter       = 100,
    converge_tol   = 1e-4,
    theta_range    = c(-4, 4),
    quad_points    = 21,
    nr_max_iter    = 20,
    nr_damp        = 1.0,
    verbose        = TRUE,
    # Priors for MAP estimation (Critical for stability)
    prior_a_mu     = 0.0,  # Log-mean (centered at a=1)
    prior_a_sd     = 0.5,  # Moderate constraint
    prior_b_mu     = 0.0,
    prior_b_sd     = 2.0,
    prior_c_alpha  = 5,    # Beta prior for guessing
    prior_c_beta   = 17,   # Mean approx 0.2
    use_priors     = TRUE
  )
  con[names(control)] <- control

  if(!is.data.frame(data)) stop("Input must be a data frame.")
  raw_data <- as.matrix(data)
  if(!is.numeric(raw_data)) stop("Input data must contain only numeric values.")

  n_persons  <- nrow(raw_data)
  n_items    <- ncol(raw_data)
  item_names <- colnames(data)
  if(is.null(item_names)) item_names <- paste0("Item_", 1:n_items)

  # Map Items to Testlets
  item_to_testlet <- rep(NA, n_items)
  for(g_idx in seq_along(group)) {
    g_items <- group[[g_idx]]
    if(is.character(g_items)) idx <- match(g_items, item_names) else idx <- g_items
    if(any(is.na(idx)) || any(idx > n_items)) stop("Invalid item indices.")
    if(any(!is.na(item_to_testlet[idx]))) stop("Items cannot overlap testlets.")
    item_to_testlet[idx] <- g_idx
  }
  n_testlets <- length(group)

  # --- USER-FRIENDLY START MESSAGE ---
  fix_discrimination <- (model == "RaschT")
  if(con$verbose) {
    cat(sprintf("\nStarting %s Estimation using %s algorithm...\n", model, method))
    if(fix_discrimination) cat("Note: Discrimination fixed to 1.0 (RaschT specification).\n")
    cat("------------------------------------------------\n")
    cat(sprintf("Data: %d Persons, %d Items, %d Testlets\n", n_persons, n_items, n_testlets))
  }

  # Helpers
  plogis_c  <- function(x) 1 / (1 + exp(-x))
  safe_prob <- function(p) {
    p[p < 1e-9] <- 1e-9;
    p[p > (1 - 1e-9)] <- 1 - 1e-9;
    p
  }

  # Initialization
  p_vals  <- colMeans(raw_data, na.rm=TRUE)
  p_vals[p_vals < 0.01] <- 0.01; p_vals[p_vals > 0.99] <- 0.99
  logit_p <- log(p_vals / (1 - p_vals))

  a_est <- rep(1.0, n_items)
  b_est <- -logit_p
  c_est <- rep(0.0, n_items)

  if(model == "RaschT") { a_est[] <- 1.0; c_est[] <- 0.0 }
  if(model == "2PLT")   { c_est[] <- 0.0 }
  if(model == "3PLT")   { c_est[] <- 0.1 }

  p_theta <- rep(0, n_persons)
  p_gamma <- matrix(0, n_persons, n_testlets)
  p_se    <- rep(NA, n_persons)
  p_gamma_se <- matrix(NA, n_persons, n_testlets)

  bad_items_idx <- c()
  is_converged  <- FALSE
  log_L_data <- NA
  iter <- 0

  # ===========================================================================
  # STRATEGY 1: EM ALGORITHM (Marginal Likelihood via Quadrature)
  # ===========================================================================
  if(method == "EM") {

    get_gh_points <- function(n) {
      if(n < 2) n <- 2
      i  <- 1:(n - 1); a <- sqrt(i / 2)
      CM <- diag(0, n); CM[cbind(1:(n-1), 2:n)] <- a; CM[cbind(2:n, 1:(n-1))] <- a
      ev <- eigen(CM, symmetric = TRUE)
      x  <- ev$values * sqrt(2); w <- (ev$vectors[1, ]^2) * sqrt(pi) / sqrt(pi)
      list(x = x[order(x)], w = w[order(x)])
    }
    gh  <- get_gh_points(con$quad_points)
    X_q <- gh$x; W_q <- gh$w; n_q <- length(X_q)

    Theta_grid <- matrix(X_q, nrow=n_q, ncol=n_q, byrow=FALSE)
    Gamma_grid <- matrix(X_q, nrow=n_q, ncol=n_q, byrow=TRUE)
    sigma_testlet <- rep(0.5, n_testlets)

    item_expected_n <- array(0, dim=c(n_items, n_q, n_q))
    item_expected_r <- array(0, dim=c(n_items, n_q, n_q))
    post_theta <- matrix(0, n_persons, n_q)

    for(it in 1:con$max_iter) {
      iter <- it

      # --- E-STEP ---
      item_expected_n[] <- 0; item_expected_r[] <- 0
      sum_E_gamma2 <- rep(0, n_testlets)
      L_testlet_theta <- array(1, dim=c(n_persons, n_testlets, n_q))

      for(t_idx in 1:n_testlets) {
        t_items <- which(item_to_testlet == t_idx)
        if(length(t_items) == 0) next
        sig_g <- sigma_testlet[t_idx]

        probs_grid <- array(0, dim=c(length(t_items), n_q, n_q))
        for(ii in seq_along(t_items)) {
          itm <- t_items[ii]
          z   <- a_est[itm] * (Theta_grid + sig_g * Gamma_grid - b_est[itm])
          probs_grid[ii,,] <- safe_prob(c_est[itm] + (1 - c_est[itm]) * plogis_c(z))
        }
        log_P <- log(probs_grid); log_Q <- log(1 - probs_grid)

        for(q_t in 1:n_q) {
          log_L_gamma <- matrix(0, nrow=n_persons, ncol=n_q)
          for(ii in seq_along(t_items)) {
            itm <- t_items[ii]; y <- raw_data[, itm]; obs <- !is.na(y)
            term <- outer(y, log_P[ii, q_t, ], "*") + outer(1-y, log_Q[ii, q_t, ], "*")
            term[!obs, ] <- 0
            log_L_gamma <- log_L_gamma + term
          }
          max_log <- apply(log_L_gamma, 1, max)
          L_gamma_scaled <- exp(log_L_gamma - max_log)
          integral <- as.vector(L_gamma_scaled %*% W_q)
          L_testlet_theta[, t_idx, q_t] <- integral * exp(max_log)
        }
      }

      L_total_theta <- apply(L_testlet_theta, c(1,3), prod)
      L_y <- as.vector(L_total_theta %*% W_q)
      log_L_data <- sum(log(L_y[L_y > 0]))

      inv_L_y <- 1/L_y; inv_L_y[L_y < 1e-300] <- 0
      post_theta <- (L_total_theta * matrix(W_q, n_persons, n_q, byrow=TRUE)) * inv_L_y

      # Accumulate Counts
      for(t_idx in 1:n_testlets) {
        t_items <- which(item_to_testlet == t_idx); if(length(t_items) == 0) next
        sig_g <- sigma_testlet[t_idx]

        probs_grid <- array(0, dim=c(length(t_items), n_q, n_q))
        for(ii in seq_along(t_items)) {
          itm <- t_items[ii]
          z   <- a_est[itm] * (Theta_grid + sig_g * Gamma_grid - b_est[itm])
          probs_grid[ii,,] <- safe_prob(c_est[itm] + (1 - c_est[itm]) * plogis_c(z))
        }
        log_P <- log(probs_grid); log_Q <- log(1 - probs_grid)

        E_g2_accum <- 0
        for(q_t in 1:n_q) {
          log_L_gamma <- matrix(0, nrow=n_persons, ncol=n_q)
          for(ii in seq_along(t_items)) {
            itm <- t_items[ii]; y <- raw_data[, itm]; obs <- !is.na(y)
            term <- outer(y, log_P[ii, q_t, ], "*") + outer(1-y, log_Q[ii, q_t, ], "*")
            term[!obs, ] <- 0
            log_L_gamma <- log_L_gamma + term
          }
          denom <- L_testlet_theta[, t_idx, q_t]; denom[denom < 1e-300] <- 1e-300
          post_g_given_t <- (exp(log_L_gamma) * matrix(W_q, n_persons, n_q, byrow=TRUE)) / denom
          joint_w <- post_g_given_t * post_theta[, q_t]
          E_g2_accum <- E_g2_accum + sum(t(joint_w) * (X_q^2))

          for(ii in seq_along(t_items)) {
            itm <- t_items[ii]; y <- raw_data[, itm]; obs <- !is.na(y)
            w_obs <- joint_w; w_obs[!obs, ] <- 0
            item_expected_n[itm, q_t, ] <- item_expected_n[itm, q_t, ] + colSums(w_obs)
            item_expected_r[itm, q_t, ] <- item_expected_r[itm, q_t, ] + colSums(w_obs * y)
          }
        }
        sum_E_gamma2[t_idx] <- E_g2_accum / n_persons
      }

      # --- M-STEP ---
      max_change <- 0
      for(i in 1:n_items) {
        if(i %in% bad_items_idx) next
        tid <- item_to_testlet[i]; sig <- if(!is.na(tid)) sigma_testlet[tid] else 0
        X_eff <- as.vector(Theta_grid + sig * Gamma_grid)
        E_n <- as.vector(item_expected_n[i,,]); E_r <- as.vector(item_expected_r[i,,])

        if(sum(E_n) < 1e-9) { bad_items_idx <- c(bad_items_idx, i); next }
        old_p <- c(a_est[i], b_est[i], c_est[i])

        for(nr in 1:con$nr_max_iter) {
          a <- a_est[i]; b <- b_est[i]; c <- c_est[i]
          lin <- a * (X_eff - b); psi <- plogis_c(lin); P <- safe_prob(c + (1-c) * psi)
          diff <- (E_r - E_n * P); w <- diff / (P * (1-P)); dP <- (1-c) * psi * (1-psi)

          ga <- sum(w * dP * (X_eff - b)); gb <- sum(w * dP * (-a)); gc <- sum(w * (1 - psi))
          w_h <- E_n / (P * (1-P))
          haa <- -sum(w_h * (dP * (X_eff - b))^2); hbb <- -sum(w_h * (dP * (-a))^2); hcc <- -sum(w_h * (1 - psi)^2)
          hab <- -sum(w_h * (dP * (X_eff - b)) * (dP * (-a))); hac <- -sum(w_h * (dP * (X_eff - b)) * (1 - psi))
          hbc <- -sum(w_h * (dP * (-a)) * (1 - psi))

          if(con$use_priors) {
            la <- log(a); pa_mu <- con$prior_a_mu; pa_sd <- con$prior_a_sd
            ga <- ga - (la - pa_mu)/(pa_sd^2 * a) - 1/a
            haa <- haa - (1 - (la - pa_mu))/(pa_sd^2 * a^2) + 1/a^2
            pb_mu <- con$prior_b_mu; pb_sd <- con$prior_b_sd
            gb <- gb - (b - pb_mu)/(pb_sd^2); hbb <- hbb - 1/(pb_sd^2)
          }

          if(model == "RaschT") { H <- matrix(hbb,1,1); G <- c(gb); idx_p <- 2 }
          else if(model %in% c("2PLT", "BiFT")) { H <- matrix(c(haa,hab,hab,hbb),2,2); G <- c(ga,gb); idx_p <- c(1,2) }
          else { H <- matrix(c(haa,hab,hac,hab,hbb,hbc,hac,hbc,hcc),3,3); G <- c(ga,gb,gc); idx_p <- 1:3 }

          diag(H) <- diag(H) - 1e-3
          delta <- tryCatch(solve(H, G), error=function(e) rep(0, length(G)))
          delta <- pmax(-0.5, pmin(0.5, delta))
          new_v <- old_p[idx_p] - con$nr_damp * delta

          if(1 %in% idx_p) new_v[which(idx_p==1)] <- max(0.1, min(4.0, new_v[which(idx_p==1)]))
          if(3 %in% idx_p) new_v[which(idx_p==3)] <- max(0.0, min(0.4, new_v[which(idx_p==3)]))

          if(model == "RaschT") b_est[i] <- new_v[1]
          else if(model %in% c("2PLT","BiFT")) { a_est[i] <- new_v[1]; b_est[i] <- new_v[2] }
          else { a_est[i] <- new_v[1]; b_est[i] <- new_v[2]; c_est[i] <- new_v[3] }
          if(max(abs(delta)) < 1e-3) break
        }
        max_change <- max(max_change, max(abs(c(a_est[i], b_est[i], c_est[i]) - old_p)))
      }

      old_sig <- sigma_testlet; sigma_testlet <- sqrt(sum_E_gamma2); sigma_testlet <- pmax(0.05, pmin(2.0, sigma_testlet))
      max_change <- max(max_change, max(abs(sigma_testlet - old_sig)))

      if(con$verbose) {
        cat(sprintf("\rIter %3d | LogLik: %10.2f | Max Change: %.6f    ", iter, log_L_data, max_change))
        flush.console()
      }

      # --- CONVERGENCE CHECK (USER MESSAGE) ---
      if(max_change < con$converge_tol) {
        if(con$verbose) cat("\n\n>>> Convergence Confirmation: Model Converged!\n")
        is_converged <- TRUE
        break
      }
    }

    # --- MAX ITER MESSAGE ---
    if(!is_converged && con$verbose) {
      cat(sprintf("\n\n>>> NOTE: Stopped because max Iteration Reached at %d.\n", con$max_iter))
      cat("    Adjust it in control=list(max_iter = ...) to have higher iteration if needed.\n")
    }

    # --- FINAL PARAMETERS MESSAGE ---
    if(con$verbose) {
      cat("\n------------------------------------------------\n")
      cat("Estimating Final Person Parameters...\n")
      cat(sprintf("Personal Parameter estimated with range [%s, %s].\n",
                  con$theta_range[1], con$theta_range[2]))
      cat("Adjust it in control=list(theta_range=c(..)) if you want a different range.\n")
    }

    p_theta <- as.vector(post_theta %*% X_q)
    p_se    <- sqrt(abs(as.vector(post_theta %*% (X_q^2)) - p_theta^2))

    # --- GAMMA POST-HOC (USER MESSAGE) ---
    if(con$verbose) cat("\nCalculating Testlet Effects (EAP) for EM...\n")

    for(t_idx in 1:n_testlets) {
      t_items <- which(item_to_testlet == t_idx); if(length(t_items) == 0) next
      sig_g <- sigma_testlet[t_idx]
      probs_grid <- array(0, dim=c(length(t_items), n_q, n_q))
      for(ii in seq_along(t_items)) {
        itm <- t_items[ii]
        z <- a_est[itm] * (Theta_grid + sig_g * Gamma_grid - b_est[itm])
        probs_grid[ii,,] <- safe_prob(c_est[itm] + (1 - c_est[itm]) * plogis_c(z))
      }
      log_P <- log(probs_grid); log_Q <- log(1 - probs_grid)
      post_gamma_agg <- matrix(0, n_persons, n_q)
      for(q_t in 1:n_q) {
        log_L <- matrix(0, n_persons, n_q)
        for(ii in seq_along(t_items)) {
          y <- raw_data[, t_items[ii]]; obs <- !is.na(y)
          term <- outer(y, log_P[ii, q_t, ], "*") + outer(1-y, log_Q[ii, q_t, ], "*")
          term[!obs, ] <- 0
          log_L <- log_L + term
        }
        denom <- L_testlet_theta[, t_idx, q_t]; denom[denom < 1e-300] <- 1e-300
        pg <- (exp(log_L) * matrix(W_q, n_persons, n_q, byrow=TRUE)) / denom
        post_gamma_agg <- post_gamma_agg + pg * post_theta[, q_t]
      }
      eg  <- as.vector(post_gamma_agg %*% X_q)
      eg2 <- as.vector(post_gamma_agg %*% (X_q^2))
      p_gamma[, t_idx] <- eg * sig_g
      p_gamma_se[, t_idx] <- sqrt(abs(eg2 - eg^2)) * sig_g
    }
  }

  # ===========================================================================
  # STRATEGY 2: MLE ALGORITHM (Joint Maximum Likelihood)
  # ===========================================================================
  else if(method == "MLE") {

    p_theta <- rowSums(raw_data, na.rm=TRUE)
    p_theta <- scale(p_theta); p_theta[is.na(p_theta)] <- 0
    p_gamma <- matrix(0, n_persons, n_testlets)

    for(it in 1:con$max_iter) {
      iter <- it

      # Update Persons
      for(p_iter in 1:2) {
        # Theta
        for(j in 1:n_persons) {
          t_val <- p_theta[j]
          eff <- rep(0, n_items)
          for(k in 1:n_testlets) {
            idx <- which(item_to_testlet == k)
            if(length(idx)>0) eff[idx] <- p_gamma[j, k]
          }
          y <- raw_data[j, ]; obs <- !is.na(y)
          if(sum(obs)==0) next
          for(nr in 1:5) {
            z <- a_est * (t_val + eff - b_est)
            p <- c_est + (1 - c_est) * plogis_c(z); p <- safe_prob(p)
            psi <- plogis_c(z); dP_dZ <- (1 - c_est) * psi * (1 - psi)
            grad <- sum((y[obs] - p[obs]) / (p[obs]*(1-p[obs])) * dP_dZ[obs] * a_est[obs])
            hess <- -sum( (dP_dZ[obs] * a_est[obs])^2 / (p[obs]*(1-p[obs])) )
            grad <- grad - t_val; hess <- hess - 1
            delta <- grad / (hess - 1e-3)
            t_val <- t_val - delta
            if(abs(delta) < 0.01) break
          }
          # Use theta_range for clamping
          p_theta[j] <- max(con$theta_range[1], min(con$theta_range[2], t_val))
        }

        # Gamma
        for(k in 1:n_testlets) {
          t_items <- which(item_to_testlet == k)
          if(length(t_items) == 0) next
          for(j in 1:n_persons) {
            g_val <- p_gamma[j, k]; t_val <- p_theta[j]
            y <- raw_data[j, t_items]; obs <- !is.na(y)
            if(sum(obs)==0) next
            for(nr in 1:5) {
              z <- a_est[t_items] * (t_val + g_val - b_est[t_items])
              psi <- plogis_c(z); p <- safe_prob(c_est[t_items] + (1 - c_est[t_items]) * psi)
              dP_dZ <- (1 - c_est[t_items]) * psi * (1 - psi)
              w <- (y[obs] - p[obs]) / (p[obs]*(1-p[obs]))
              grad <- sum(w * dP_dZ[obs] * a_est[t_items])
              hess <- -sum( (dP_dZ[obs] * a_est[t_items])^2 / (p[obs]*(1-p[obs])) )
              grad <- grad - g_val/0.25; hess <- hess - 1/0.25
              delta <- grad / (hess - 1e-3)
              g_val <- g_val - delta
              if(abs(delta) < 0.01) break
            }
            p_gamma[j, k] <- max(-2, min(2, g_val))
          }
        }
      }
      p_theta <- as.vector(scale(p_theta))

      # Update Items
      max_change <- 0
      for(i in 1:n_items) {
        if(i %in% bad_items_idx) next
        tid <- item_to_testlet[i]; eff <- if(!is.na(tid)) p_gamma[, tid] else rep(0, n_persons)
        X_eff <- p_theta + eff
        y_vec <- raw_data[, i]; obs <- !is.na(y_vec)
        old_p <- c(a_est[i], b_est[i], c_est[i])

        for(nr in 1:con$nr_max_iter) {
          a <- a_est[i]; b <- b_est[i]; c <- c_est[i]
          X_o <- X_eff[obs]; y_o <- y_vec[obs]
          lin <- a * (X_o - b); psi <- plogis_c(lin); P <- safe_prob(c + (1-c) * psi)
          diff <- (y_o - P); w <- diff / (P * (1-P)); dP <- (1-c) * psi * (1-psi)
          ga <- sum(w * dP * (X_o - b)); gb <- sum(w * dP * (-a)); gc <- sum(w * (1 - psi))
          w_h <- 1 / (P * (1-P))
          haa <- -sum(w_h * (dP * (X_o - b))^2); hbb <- -sum(w_h * (dP * (-a))^2); hcc <- -sum(w_h * (1 - psi)^2)
          hab <- -sum(w_h * (dP * (X_o - b)) * (dP * (-a))); hac <- -sum(w_h * (dP * (X_o - b)) * (1 - psi))
          hbc <- -sum(w_h * (dP * (-a)) * (1 - psi))

          if(con$use_priors) {
            la <- log(a); pa_mu <- con$prior_a_mu; pa_sd <- con$prior_a_sd
            ga <- ga - (la - pa_mu)/(pa_sd^2 * a) - 1/a
            haa <- haa - (1 - (la - pa_mu))/(pa_sd^2 * a^2) + 1/a^2
            pb_mu <- con$prior_b_mu; pb_sd <- con$prior_b_sd
            gb <- gb - (b - pb_mu)/(pb_sd^2); hbb <- hbb - 1/(pb_sd^2)
          }

          if(model == "RaschT") { H <- matrix(hbb,1,1); G <- c(gb); idx_p <- 2 }
          else if(model %in% c("2PLT", "BiFT")) { H <- matrix(c(haa,hab,hab,hbb),2,2); G <- c(ga,gb); idx_p <- c(1,2) }
          else { H <- matrix(c(haa,hab,hac,hab,hbb,hbc,hac,hbc,hcc),3,3); G <- c(ga,gb,gc); idx_p <- 1:3 }

          diag(H) <- diag(H) - 1e-3
          delta <- tryCatch(solve(H, G), error=function(e) rep(0, length(G)))
          delta <- pmax(-0.5, pmin(0.5, delta))

          new_v <- old_p[idx_p] - con$nr_damp * delta
          if(1 %in% idx_p) new_v[which(idx_p==1)] <- max(0.1, min(4.0, new_v[which(idx_p==1)]))
          if(3 %in% idx_p) new_v[which(idx_p==3)] <- max(0.0, min(0.4, new_v[which(idx_p==3)]))

          if(model == "RaschT") b_est[i] <- new_v[1]
          else if(model %in% c("2PLT","BiFT")) { a_est[i] <- new_v[1]; b_est[i] <- new_v[2] }
          else { a_est[i] <- new_v[1]; b_est[i] <- new_v[2]; c_est[i] <- new_v[3] }
          if(max(abs(delta)) < 1e-3) break
        }
        max_change <- max(max_change, max(abs(c(a_est[i], b_est[i], c_est[i]) - old_p)))
      }

      curr_ll <- 0
      for(i in 1:n_items) {
        tid <- item_to_testlet[i]; eff <- if(!is.na(tid)) p_gamma[, tid] else rep(0, n_persons)
        z <- a_est[i] * (p_theta + eff - b_est[i])
        p <- safe_prob(c_est[i] + (1 - c_est[i]) * plogis_c(z))
        y <- raw_data[, i]; obs <- !is.na(y)
        curr_ll <- curr_ll + sum(y[obs]*log(p[obs]) + (1-y[obs])*log(1-p[obs]))
      }
      log_L_data <- curr_ll

      if(con$verbose) {
        cat(sprintf("\rIter %3d | LogLik: %10.2f | Max Change: %.6f    ", iter, log_L_data, max_change))
        flush.console()
      }

      # --- CONVERGENCE CHECK (USER MESSAGE) ---
      if(max_change < con$converge_tol) {
        if(con$verbose) cat("\n\n>>> Convergence Confirmation: Model Converged!\n")
        is_converged <- TRUE
        break
      }
    }

    # --- MAX ITER MESSAGE ---
    if(!is_converged && con$verbose) {
      cat(sprintf("\n\n>>> NOTE: Stopped because max Iteration Reached at %d.\n", con$max_iter))
      cat("    Adjust it in control=list(max_iter = ...) to have higher iteration if needed.\n")
    }
    if(con$verbose) {
      cat("\n------------------------------------------------\n")
      cat("Estimating Final Person Parameters...\n")
      cat(sprintf("Personal Parameter estimated with range [%s, %s].\n",
                  con$theta_range[1], con$theta_range[2]))
      cat("Adjust it in control=list(theta_range=c(..)) if you want a different range.\n")
    }

    # Approx SEs for MLE
    for(j in 1:n_persons) {
      eff <- rep(0, n_items)
      for(k in 1:n_testlets) if(!is.na(match(k, item_to_testlet))) eff[item_to_testlet==k] <- p_gamma[j,k]
      z <- a_est * (p_theta[j] + eff - b_est)
      psi <- plogis_c(z); dP_dZ <- (1 - c_est) * psi * (1 - psi)
      p <- safe_prob(c_est + (1 - c_est) * psi)
      info <- sum((dP_dZ * a_est)^2 / (p*(1-p)), na.rm=TRUE)
      p_se[j] <- 1/sqrt(info + 1)
    }
    for(k in 1:n_testlets) {
      t_items <- which(item_to_testlet == k)
      if(length(t_items)==0) next
      for(j in 1:n_persons) {
        z <- a_est[t_items] * (p_theta[j] + p_gamma[j,k] - b_est[t_items])
        psi <- plogis_c(z); dP_dZ <- (1 - c_est[t_items]) * psi * (1 - psi)
        p <- safe_prob(c_est[t_items] + (1 - c_est[t_items]) * psi)
        info <- sum((dP_dZ * a_est[t_items])^2 / (p*(1-p)), na.rm=TRUE)
        p_gamma_se[j, k] <- 1/sqrt(info + 4)
      }
    }
  }

  # ===========================================================================
  # FINAL SEs (ITEMS) & OUTPUT
  # ===========================================================================
  se_a <- rep(NA, n_items); se_b <- rep(NA, n_items); se_c <- rep(NA, n_items)

  for(i in 1:n_items) {
    if(i %in% bad_items_idx) next

    tid <- item_to_testlet[i];
    if(method=="EM") {
      sig <- if(!is.na(tid)) sigma_testlet[tid] else 0
      X_eff <- as.vector(Theta_grid + sig * Gamma_grid)
      E_n <- as.vector(item_expected_n[i,,])
      weights <- E_n
    } else {
      eff <- if(!is.na(tid)) p_gamma[, tid] else 0
      X_eff <- (p_theta + eff)[!is.na(raw_data[,i])]
      weights <- 1
    }

    a <- a_est[i]; b <- b_est[i]; c <- c_est[i]
    lin <- a * (X_eff - b); psi <- plogis_c(lin); P <- safe_prob(c + (1-c) * psi)
    dP <- (1-c) * psi * (1-psi)
    w_h <- weights / (P * (1-P))

    haa <- -sum(w_h * (dP * (X_eff - b))^2); hbb <- -sum(w_h * (dP * (-a))^2)
    hcc <- -sum(w_h * (1 - psi)^2); hab <- -sum(w_h * (dP * (X_eff - b)) * (dP * (-a)))
    hac <- -sum(w_h * (dP * (X_eff - b)) * (1 - psi)); hbc <- -sum(w_h * (dP * (-a)) * (1 - psi))

    if(con$use_priors) {
      la <- log(a); pa_mu <- con$prior_a_mu; pa_sd <- con$prior_a_sd
      haa <- haa - (1 - (la - pa_mu))/(pa_sd^2 * a^2) + 1/a^2
      pb_mu <- con$prior_b_mu; pb_sd <- con$prior_b_sd
      hbb <- hbb - 1/(pb_sd^2)
    }

    H <- NULL; idx_map <- NULL
    if(model == "RaschT") { H <- matrix(hbb,1,1); idx_map <- 2 }
    else if(model %in% c("2PLT", "BiFT")) { H <- matrix(c(haa,hab,hab,hbb),2,2); idx_map <- c(1,2) }
    else { H <- matrix(c(haa,hab,hac,hab,hbb,hbc,hac,hbc,hcc),3,3); idx_map <- 1:3 }

    cv <- tryCatch(solve(-H), error=function(e) matrix(NA, nrow(H), ncol(H)))
    s <- sqrt(diag(cv))
    if(1 %in% idx_map) se_a[i] <- s[which(idx_map==1)]
    if(2 %in% idx_map) se_b[i] <- s[which(idx_map==2)]
    if(3 %in% idx_map) se_c[i] <- s[which(idx_map==3)]
  }

  if(model == "RaschT") { se_a[] <- 0.0; se_c[] <- 0.0 }
  else if(model %in% c("2PLT", "BiFT")) { se_c[] <- 0.0 }

  out_items <- data.frame(
    item = item_names,
    discrimination = a_est, discrimination_se = se_a,
    difficulty = b_est, difficulty_se = se_b,
    guessing = c_est, guessing_se = se_c
  )
  if(length(bad_items_idx) > 0) out_items[bad_items_idx, -1] <- NA
  rownames(out_items) <- NULL

  out_persons <- data.frame(person=1:n_persons, ability=p_theta, ability_se=p_se)

  t_cols <- data.frame(matrix(NA, nrow=n_persons, ncol=2*n_testlets))
  t_names <- c()
  for(k in 1:n_testlets) {
    t_cols[, (k*2)-1] <- p_gamma[, k]
    t_cols[, (k*2)]   <- p_gamma_se[, k]
    t_names <- c(t_names, paste0("testlet_", k), paste0("testlet_", k, "_se"))
  }
  colnames(t_cols) <- t_names
  out_persons <- cbind(out_persons, t_cols)
  rownames(out_persons) <- NULL

  n_p <- switch(model, "RaschT"=n_items+n_testlets, "2PLT"=2*n_items+n_testlets, "BiFT"=2*n_items+n_testlets, "3PLT"=3*n_items+n_testlets)
  out_fit <- data.frame(LogLikelihood=log_L_data, AIC=2*n_p - 2*log_L_data, BIC=n_p*log(n_persons) - 2*log_L_data, Iterations=iter)
  rownames(out_fit) <- NULL

  se_cols_items <- grep("_se$", names(out_items))
  se_cols_persons <- grep("_se$", names(out_persons))

  na_items <- any(is.na(out_items[, se_cols_items]))
  na_persons <- any(is.na(out_persons[, se_cols_persons]))

  if(na_items || na_persons) {
    parts <- c()
    if(na_items) parts <- c(parts, "Item Parameters")
    if(na_persons) parts <- c(parts, "Person Parameters")
    msg <- paste0("\nWarning: Standard error detected NA in ", paste(parts, collapse=" & "),
                  ", indicating estimation for this is not stable (Hessian Singular).")
    message(msg)
  }

  if(con$verbose) {
    cat("Finished All Estimation.\n")
    cat("------------------------------------------------\n")
  }

  list(item_params=out_items, person_params=out_persons, model_fit=out_fit)
}
