#' Mixed Item Response Model Estimation (Dichotomous & Polytomous)
#'
#' @description
#' Provides a estimation framework for a broad class of different item response theory models.
#' This function can model different combinations of item categories.
#'
#' @param data A N x J data.frame. Binary items must be 0/1. Polytomous items should be continuous integers (0, 1, 2...).
#' @param model A character vector of length J (one model per item).
#'              Supported: "Rasch", "2PL" (2-Parameter Logistic), "3PL" (3-Parameter Logistic), "GRM" (Graded Response Model), "GPCM" (Generalized Partial Credit Model), "PCM" (Partial Credit Model).
#'              If a single string is provided, it is applied to all items.
#' @param method String. "EM" (Marginal Maximum Likelihood via Expectation-Maximization) or "MLE" (Joint Maximum Likelihood).
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
#'   \item \code{item_params}: Data frame of item parameters (discrimination, difficulty/thresholds, guessing).
#'   \item \code{person_params}: A data frame of estimated person abilities (theta) and standard errors.
#'   \item \code{model_fit}: A data frame containing fit statistics such as Akaike’s Information Criterion (AIC) and the Bayesian Information Criterion (BIC).
#'   \item \code{settings}: A list of control parameters used in the estimation.
#' }
#' @examples
#'   # --- Example 1: Simulation (Mixed 2PL + GPCM) ---
#'   set.seed(2025)
#'   N <- 100
#'   n_bin <- 5
#'   n_poly <- 2
#'   J <- n_bin + n_poly
#'
#'   # 1. Generate Theta (Wide range to match user request)
#'   true_theta <- rnorm(N, mean = 0, sd = 3)
#'
#'   # 2. Simulation Helper: GPCM
#'   sim_gpcm <- function(theta, a, steps) {
#'     n_cat <- length(steps) + 1
#'     probs <- matrix(0, length(theta), n_cat)
#'     for(k in 1:n_cat) {
#'       score <- k - 1
#'       if(score == 0) numer <- rep(0, length(theta))
#'       else numer <- a * (score * theta - sum(steps[1:score]))
#'       probs[, k] <- exp(numer)
#'     }
#'     probs <- probs / rowSums(probs)
#'     apply(probs, 1, function(p) sample(0:(n_cat-1), 1, prob=p))
#'   }
#'
#'   # 3. Create Data
#'   data_sim <- data.frame(matrix(NA, nrow = N, ncol = J))
#'   colnames(data_sim) <- paste0("Item_", 1:J)
#'
#'   # Binary Items (2PL)
#'   a_bin <- runif(n_bin, 0.8, 1.5)
#'   b_bin <- seq(-3, 3, length.out = n_bin)
#'   for(j in 1:n_bin) {
#'     prob <- 1 / (1 + exp(-(a_bin[j] * (true_theta - b_bin[j]))))
#'     data_sim[, j] <- rbinom(N, 1, prob)
#'   }
#'
#'   # Polytomous Items (GPCM)
#'   # Item 6: 2 steps (-2, 2)
#'   data_sim[, 6] <- sim_gpcm(true_theta, a=1.0, steps=c(-2, 2))
#'   # Item 7: 5 steps
#'   data_sim[, 7] <- sim_gpcm(true_theta, a=1.2, steps=c(-5, -2.5, 0, 2.5, 5))
#'
#'   # 4. Run Estimation
#'   # Note: Wide theta_range needed due to SD=3 in simulation
#'   my_models <- c(rep("2PL", n_bin), rep("GPCM", n_poly))
#'
#'   res <- mixed_irt(data = data_sim, model = my_models, method = "EM",
#'                    control = list(max_iter = 20, theta_range = c(-6, 6)))
#'
#'   head(res$item_params)
#'   print(res$model_fit)
#'   \donttest{
#'   # --- Example 2: With Package Data ---
#'   data("ela2", package = "tirt")
#'
#'   # Define Models (7 Binary, 3 Poly)
#'   real_models <- c(rep("2PL", 7), rep("GRM", 3))
#'
#'   # Run Estimation
#'   real_res <- mixed_irt(ela2, model = real_models, method = "EM",
#'                         control = list(max_iter = 10))
#'
#'   head(real_res$item_params)
#'   print(real_res$model_fit)
#'   }
#' @export
mixed_irt <- function(data,
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
    nr_damp = 0.8,
    verbose = TRUE
  )
  con[names(control)] <- control

  # --- Input Validation ---
  if(!is.data.frame(data)) stop("Input must be a data frame.")
  raw_mat <- as.matrix(data)
  if(!is.numeric(raw_mat)) stop("Input data must be numeric.")

  N <- nrow(raw_mat)
  J <- ncol(raw_mat)
  item_names <- colnames(data)
  if(is.null(item_names)) item_names <- paste0("Item_", 1:J)

  # Expand model if single string
  if(length(model) == 1) model <- rep(model, J)
  if(length(model) != J) stop("Length of 'model' vector must match number of items.")

  # Define Model Types
  binary_models <- c("Rasch", "2PL", "3PL")
  poly_models <- c("GRM", "GPCM", "PCM")

  # --- Data Preprocessing ---
  n_cats <- numeric(J)
  for(j in 1:J) {
    if(model[j] %in% poly_models) {
      min_val <- min(raw_mat[,j], na.rm=TRUE)
      if(is.finite(min_val) && min_val > 0) {
        raw_mat[,j] <- raw_mat[,j] - min_val
      }
      n_cats[j] <- max(raw_mat[,j], na.rm=TRUE) + 1
    } else {
      n_cats[j] <- 2
    }
  }

  max_K <- max(n_cats)

  # Identify valid items (Variance > 0)
  item_vars <- apply(raw_mat, 2, var, na.rm=TRUE)
  valid_items_idx <- which(item_vars > 0)
  bad_items_idx <- setdiff(1:J, valid_items_idx)

  # --- Initialize Parameters ---
  a <- rep(1, J)
  b <- rep(0, J)
  g <- rep(0, J)
  d <- matrix(0, J, max_K)

  # Initialize Standard Errors (SE)
  se_a <- rep(NA, J)
  se_b <- rep(NA, J)
  se_g <- rep(NA, J)
  se_d <- matrix(NA, J, max_K)

  for(j in valid_items_idx) {
    valid_resp <- raw_mat[!is.na(raw_mat[,j]), j]

    if(model[j] %in% binary_models) {
      p_val <- mean(valid_resp)
      p_bounded <- pmin(pmax(p_val, 0.01), 0.99)
      b[j] <- -log(p_bounded / (1 - p_bounded))
      if(model[j] == "3PL") g[j] <- 0.2
    } else {
      props <- table(factor(valid_resp, levels=0:(n_cats[j]-1))) / length(valid_resp)
      cum_props <- cumsum(props)
      cum_props <- pmin(pmax(cum_props, 0.01), 0.99)
      thresh <- -qlogis(cum_props[-length(cum_props)])

      if(model[j] == "GRM") {
        d[j, 1:(n_cats[j]-1)] <- sort(thresh)
      } else {
        K <- n_cats[j]
        steps <- seq(-0.5, 0.5, length.out = K-1)
        d[j, 2:K] <- steps
      }
    }
  }

  # Initialize Theta
  raw_scores <- rowMeans(raw_mat, na.rm=TRUE)
  theta <- as.vector(scale(raw_scores))
  theta[is.na(theta)] <- 0

  # Quadrature Setup
  nodes <- seq(con$theta_range[1], con$theta_range[2], length.out = con$quad_points)
  weights <- dnorm(nodes)
  weights <- weights / sum(weights)

  # --- 2. Helper Functions ---

  get_probs_j <- function(th_vec, j_idx) {
    mod <- model[j_idx]
    if(mod %in% binary_models) {
      z <- a[j_idx] * (th_vec - b[j_idx])
      z <- pmin(pmax(z, -30), 30)
      p <- g[j_idx] + (1 - g[j_idx]) * (1 / (1 + exp(-z)))
      return(cbind(1-p, p))
    } else {
      n_cat <- n_cats[j_idx]
      n_q <- length(th_vec)
      probs <- matrix(0, n_q, n_cat)
      if (mod %in% c("GPCM", "PCM")) {
        for(k in 1:n_cat) {
          z <- a[j_idx] * (k-1) * th_vec - d[j_idx, k]
          z <- pmin(pmax(z, -50), 50)
          probs[,k] <- exp(z)
        }
        probs <- probs / rowSums(probs)
      } else { # GRM
        P_star <- matrix(0, n_q, n_cat + 1)
        P_star[,1] <- 1; P_star[,n_cat + 1] <- 0
        eff_d <- d[j_idx, 1:(n_cat-1)]
        if(is.unsorted(eff_d)) eff_d <- sort(eff_d)
        for(k in 1:(n_cat-1)) {
          z <- a[j_idx] * (th_vec - eff_d[k])
          z <- pmin(pmax(z, -30), 30)
          P_star[, k+1] <- 1 / (1 + exp(-z))
        }
        for(k in 1:n_cat) probs[,k] <- P_star[,k] - P_star[,k+1]
        probs <- pmax(probs, 1e-10)
      }
      return(probs)
    }
  }

  optim_binary <- function(r_vec, n_vec, th_vec, cur_a, cur_b, cur_g, mod_type) {
    get_P_bin <- function(t, pa, pb, pg) {
      z <- pa * (t - pb); z <- pmin(pmax(z, -30), 30)
      pg + (1 - pg) * (1 / (1 + exp(-z)))
    }
    if(mod_type == "Rasch") active <- c(FALSE, TRUE, FALSE)
    else if(mod_type == "2PL") active <- c(TRUE, TRUE, FALSE)
    else active <- c(TRUE, TRUE, TRUE)

    par_vec <- c(cur_a, cur_b, cur_g)[active]

    ll_func <- function(p_v) {
      full_p <- c(1, 0, 0); full_p[active] <- p_v
      if(!active[1]) full_p[1] <- 1
      P_tmp <- get_P_bin(th_vec, full_p[1], full_p[2], full_p[3])
      P_tmp <- pmin(pmax(P_tmp, 1e-9), 1-1e-9)
      sum(r_vec * log(P_tmp) + (n_vec - r_vec) * log(1 - P_tmp))
    }

    for(k in 1:con$nr_max_iter) {
      f0 <- ll_func(par_vec); h <- 1e-4
      grad <- numeric(length(par_vec)); hess <- matrix(0, length(par_vec), length(par_vec))
      for(i in 1:length(par_vec)) {
        tmp <- par_vec; tmp[i] <- tmp[i] + h; grad[i] <- (ll_func(tmp) - f0)/h
      }
      for(i in 1:length(par_vec)) {
        tmp_kk <- par_vec; tmp_kk[i] <- tmp_kk[i] + h; tmp_mk <- par_vec; tmp_mk[i] <- tmp_mk[i] - h
        hess[i,i] <- (ll_func(tmp_kk) - 2*f0 + ll_func(tmp_mk)) / h^2
      }
      diag(hess) <- diag(hess) - 1e-5
      step <- try(solve(hess, grad), silent=TRUE)
      if(inherits(step, "try-error")) break
      par_vec <- par_vec - step * con$nr_damp
      if(active[1]) par_vec[1] <- max(par_vec[1], 0.01)
      if(active[3]) par_vec[3] <- min(max(par_vec[3], 0), 0.4)
      if(max(abs(step)) < 1e-5) break
    }

    se_vec <- tryCatch(sqrt(diag(solve(-hess))), error=function(e) rep(NA, length(par_vec)))

    res <- list(a=if(active[1]) par_vec[1] else 1, b=0, g=0, se_a=NA, se_b=NA, se_g=NA)
    idx <- 1
    if(active[1]) { res$se_a <- se_vec[idx]; idx <- idx+1 }
    if(active[2]) { res$b <- par_vec[idx]; res$se_b <- se_vec[idx]; idx <- idx+1 }
    if(active[3]) { res$g <- par_vec[idx]; res$se_g <- se_vec[idx] }
    return(res)
  }

  optim_poly <- function(r_mat, n_vec, th_vec, cur_a, cur_d, mod_type, n_cat) {
    if(mod_type == "PCM") param_vec <- cur_d[2:n_cat]
    else if(mod_type == "GRM") param_vec <- c(cur_a, cur_d[1:(n_cat-1)])
    else param_vec <- c(cur_a, cur_d[2:n_cat])

    ll_func <- function(p_v) {
      l_a <- if(mod_type=="PCM") 1 else p_v[1]
      l_d <- numeric(max_K)
      if(mod_type=="PCM") l_d[2:n_cat] <- p_v
      else if(mod_type=="GRM") l_d[1:(n_cat-1)] <- p_v[-1]
      else l_d[2:n_cat] <- p_v[-1]

      n_q <- length(th_vec); probs <- matrix(0, n_q, n_cat)
      if (mod_type %in% c("GPCM", "PCM")) {
        for(k in 1:n_cat) {
          z <- l_a * (k-1) * th_vec - l_d[k]; z <- pmin(pmax(z, -50), 50); probs[,k] <- exp(z)
        }
        probs <- probs / rowSums(probs)
      } else {
        P_s <- matrix(0, n_q, n_cat+1); P_s[,1] <- 1; ed <- l_d[1:(n_cat-1)]
        if(is.unsorted(ed)) ed <- sort(ed)
        for(k in 1:(n_cat-1)) {
          z <- l_a * (th_vec - ed[k]); z <- pmin(pmax(z, -30), 30); P_s[, k+1] <- 1/(1+exp(-z))
        }
        for(k in 1:n_cat) probs[,k] <- P_s[,k] - P_s[,k+1]
        probs <- pmax(probs, 1e-10)
      }
      sum(r_mat * log(probs))
    }

    for(k in 1:con$nr_max_iter) {
      f0 <- ll_func(param_vec); h <- 1e-4
      grad <- numeric(length(param_vec)); hess <- matrix(0, length(param_vec), length(param_vec))
      for(i in 1:length(param_vec)) {
        t <- param_vec; t[i] <- t[i]+h; grad[i] <- (ll_func(t)-f0)/h
      }
      for(i in 1:length(param_vec)) {
        tp <- param_vec; tp[i] <- tp[i]+h; tm <- param_vec; tm[i] <- tm[i]-h
        hess[i,i] <- (ll_func(tp)-2*f0+ll_func(tm))/h^2
      }
      diag(hess) <- diag(hess) - 1e-4
      step <- try(solve(hess, grad), silent=TRUE)
      if(inherits(step, "try-error")) break
      param_vec <- param_vec - step * con$nr_damp
      if(mod_type != "PCM") param_vec[1] <- max(param_vec[1], 0.05)
      if(max(abs(step)) < 1e-4) break
    }

    se_vec <- tryCatch(sqrt(diag(solve(-hess))), error=function(e) rep(NA, length(param_vec)))

    ret_a <- if(mod_type=="PCM") 1 else param_vec[1]; ret_d <- cur_d
    ret_se_a <- NA; ret_se_d <- rep(NA, max_K)

    if(mod_type=="PCM") {
      ret_d[2:n_cat] <- param_vec; ret_se_d[2:n_cat] <- se_vec
    } else if(mod_type=="GRM") {
      ret_d[1:(n_cat-1)] <- sort(param_vec[-1]); ret_se_a <- se_vec[1]; ret_se_d[1:(n_cat-1)] <- se_vec[-1]
    } else {
      ret_d[2:n_cat] <- param_vec[-1]; ret_se_a <- se_vec[1]; ret_se_d[2:n_cat] <- se_vec[-1]
    }
    return(list(a=ret_a, d=ret_d, se_a=ret_se_a, se_d=ret_se_d))
  }

  # --- 3. Main Estimation Loop ---

  if(con$verbose) cat(sprintf("\nStarting Mixed IRT (%s) Estimation...\n", method))
  is_converged <- FALSE

  for(iter in 1:con$max_iter) {

    # E-Step / MLE Person Step
    if(method == "EM") {
      L_iq <- matrix(0, N, con$quad_points)
      for(j in valid_items_idx) {
        probs_jq <- get_probs_j(nodes, j)
        resp <- raw_mat[, j]; valid_rows <- !is.na(resp)
        if(model[j] %in% binary_models) {
          p1 <- pmax(probs_jq[,2], 1e-10); p0 <- pmax(probs_jq[,1], 1e-10)
          is_1 <- which(valid_rows & resp == 1); is_0 <- which(valid_rows & resp == 0)
          if(length(is_1)>0) L_iq[is_1,] <- L_iq[is_1,] + matrix(log(p1), length(is_1), con$quad_points, byrow=TRUE)
          if(length(is_0)>0) L_iq[is_0,] <- L_iq[is_0,] + matrix(log(p0), length(is_0), con$quad_points, byrow=TRUE)
        } else {
          p_indices <- resp[valid_rows] + 1
          for(k in 1:n_cats[j]) {
            is_cat <- which(p_indices == k)
            if(length(is_cat)>0) L_iq[valid_rows,][is_cat,] <- L_iq[valid_rows,][is_cat,] + matrix(log(pmax(probs_jq[,k], 1e-10)), length(is_cat), con$quad_points, byrow=TRUE)
          }
        }
      }
      F_iq <- exp(L_iq) * matrix(weights, N, con$quad_points, byrow=TRUE)
      Posterior_iq <- F_iq / rowSums(F_iq); Posterior_iq[is.na(Posterior_iq)] <- 0
    } else {
      # MLE Theta update
      for(i in 1:N) {
        valid <- !is.na(raw_mat[i,]) & (1:J %in% valid_items_idx); if(sum(valid)==0) next
        ti <- theta[i]
        for(nr in 1:5) {
          num <- 0; den <- 0
          for(j in which(valid)) {
            x <- raw_mat[i,j]; probs <- get_probs_j(ti, j)
            if(model[j] %in% binary_models) {
              p <- probs[1,2]; num <- num + a[j]*(x - p); den <- den - (a[j]^2)*p*(1-p)
            } else {
              h <- 1e-4; p0 <- max(probs[1, x+1], 1e-10)
              pp <- get_probs_j(ti+h, j)[1, x+1]; pm <- get_probs_j(ti-h, j)[1, x+1]
              d1 <- (log(pp)-log(pm))/(2*h); d2 <- (log(pp)-2*log(p0)+log(pm))/(h^2)
              num <- num + d1; den <- den + d2
            }
          }
          if(abs(den) < 1e-6) break; step <- num/den; ti <- ti - step
          ti <- max(min(ti, con$theta_range[2]), con$theta_range[1]); if(abs(step) < 1e-3) break
        }
        theta[i] <- ti
      }
      theta <- theta - mean(theta)
    }

    # M-Step / MLE Item Update
    max_change <- 0
    for(j in valid_items_idx) {
      old_par <- c(a[j], b[j], g[j], d[j,])
      if(method == "EM") {
        use_nodes <- nodes
        if(model[j] %in% binary_models) {
          valid_rows <- !is.na(raw_mat[,j])
          r_jq <- colSums(Posterior_iq[valid_rows,] * raw_mat[valid_rows, j])
          n_jq <- colSums(Posterior_iq[valid_rows,])
          res <- optim_binary(r_jq, n_jq, use_nodes, a[j], b[j], g[j], model[j])
          a[j] <- res$a; b[j] <- res$b; g[j] <- res$g
          se_a[j] <- res$se_a; se_b[j] <- res$se_b; se_g[j] <- res$se_g
        } else {
          R_qk <- matrix(0, con$quad_points, n_cats[j])
          valid_rows <- !is.na(raw_mat[,j]); resp <- raw_mat[valid_rows, j]
          for(k in 1:n_cats[j]) {
            idx <- which(resp == (k-1))
            if(length(idx)>0) { full_idx <- which(valid_rows)[idx]; R_qk[,k] <- colSums(Posterior_iq[full_idx, , drop=FALSE]) }
          }
          res <- optim_poly(R_qk, rowSums(R_qk), use_nodes, a[j], d[j,], model[j], n_cats[j])
          a[j] <- res$a; d[j,] <- res$d
          se_a[j] <- res$se_a; se_d[j,] <- res$se_d
        }
      } else {
        valid <- !is.na(raw_mat[,j]); use_nodes <- theta[valid]; y <- raw_mat[valid, j]
        if(model[j] %in% binary_models) {
          res <- optim_binary(y, rep(1, sum(valid)), use_nodes, a[j], b[j], g[j], model[j])
          a[j] <- res$a; b[j] <- res$b; g[j] <- res$g
          se_a[j] <- res$se_a; se_b[j] <- res$se_b; se_g[j] <- res$se_g
        } else {
          n_resp <- length(y); R_qk <- matrix(0, n_resp, n_cats[j]); R_qk[cbind(1:n_resp, y+1)] <- 1
          res <- optim_poly(R_qk, rep(1, n_resp), use_nodes, a[j], d[j,], model[j], n_cats[j])
          a[j] <- res$a; d[j,] <- res$d
          se_a[j] <- res$se_a; se_d[j,] <- res$se_d
        }
      }
      max_change <- max(max_change, abs(c(a[j], b[j], g[j], d[j,]) - old_par), na.rm=TRUE)
    }

    if(con$verbose) cat(paste0("\rIteration ", iter, ": Max Param Change = ", round(max_change, 5), "   "))

    if(max_change < con$converge_tol) {
      if(con$verbose) cat("\n\n>>> Convergence Confirmation: Model Converged!\n")
      is_converged <- TRUE
      break
    }
  }

  if(!is_converged && con$verbose) {
    cat(sprintf("\n\n>>> NOTE: Stopped because max Iteration Reached at %d.\n", con$max_iter))
    cat("    Adjust it in control=list(max_iter = ...) to have higher iteration if needed.\n")
  }

  # --- 4. Final Person Estimation & Formatting ---

  if(con$verbose) {
    cat("\n------------------------------------------------\n")
    cat("Estimating Final Person Parameters...\n")
    cat(sprintf("Personal Parameter estimated with range [%s, %s].\n", con$theta_range[1], con$theta_range[2]))
    cat("Adjust it in control=list(theta_range=c(..)) if you want a different range.\n")
  }

  final_theta <- numeric(N)
  final_se <- numeric(N)
  log_lik <- 0

  for(i in 1:N) {
    valid <- !is.na(raw_mat[i,]) & (1:J %in% valid_items_idx)
    if(sum(valid)==0) { final_theta[i] <- NA; next }

    ti <- if(method=="EM") 0 else theta[i]
    for(k in 1:20) {
      num <- 0; den <- 0
      for(j in which(valid)) {
        x <- raw_mat[i,j]; probs <- get_probs_j(ti, j)
        if(model[j] %in% binary_models) {
          p <- probs[1,2]; num <- num + a[j]*(x - p); den <- den - (a[j]^2)*p*(1-p)
          if(k==20) log_lik <- log_lik + (x*log(p) + (1-x)*log(1-p))
        } else {
          h <- 1e-4; p0 <- max(probs[1, x+1], 1e-10)
          pp <- get_probs_j(ti+h, j)[1, x+1]; pm <- get_probs_j(ti-h, j)[1, x+1]
          d1 <- (log(pp)-log(pm))/(2*h); d2 <- (log(pp)-2*log(p0)+log(pm))/(h^2)
          num <- num + d1; den <- den + d2
          if(k==20) log_lik <- log_lik + log(p0)
        }
      }
      if(abs(den) < 1e-6) break; step <- num/den; step <- max(min(step, 1.0), -1.0)
      ti <- ti - step; if(abs(step) < 1e-4) break
    }
    final_theta[i] <- max(min(ti, con$theta_range[2]), con$theta_range[1])
    if(!is.na(den) && abs(den) > 1e-6) final_se[i] <- 1/sqrt(abs(den)) else final_se[i] <- NA
  }

  if(con$verbose) cat("Person Parameter Estimation Finished.\n")

  # --- 5. Output Construction with Descriptives & SEs ---

  # Initial Data Frame
  out_items <- data.frame(item = item_names,
                          model = model,
                          stringsAsFactors = FALSE)

  # Interleave Estimates and SEs
  out_items$discrimination <- round(a, 3)
  out_items$discrimination_se <- round(se_a, 3)

  out_items$difficulty <- ifelse(model %in% binary_models, round(b, 3), NA)
  out_items$difficulty_se <- ifelse(model %in% binary_models, round(se_b, 3), NA)

  out_items$guess <- ifelse(model == "3PL", round(g, 3), NA)
  out_items$guess_se <- ifelse(model == "3PL", round(se_g, 3), NA)

  for(k in 1:(max_K-1)) {
    col_name <- paste0("threshold_", k)
    col_name_se <- paste0("threshold_", k, "_se")

    vals <- numeric(J)
    vals_se <- numeric(J)

    for(j in 1:J) {
      if(model[j] %in% poly_models) {
        if(model[j] == "GRM") {
          vals[j] <- d[j, k]
          vals_se[j] <- se_d[j, k]
        } else {
          if((k+1) <= n_cats[j]) {
            vals[j] <- d[j, k+1]
            vals_se[j] <- se_d[j, k+1]
          } else {
            vals[j] <- NA; vals_se[j] <- NA
          }
        }
      } else {
        vals[j] <- NA; vals_se[j] <- NA
      }
    }
    out_items[[col_name]] <- round(vals, 3)
    out_items[[col_name_se]] <- round(vals_se, 3)
  }

  # Calculate Descriptive Statistics and append as LAST columns
  n_count <- colSums(!is.na(raw_mat))
  p_val <- colMeans(raw_mat, na.rm = TRUE)

  out_items$number <- n_count
  out_items$pvalue <- round(p_val, 3)

  if(length(bad_items_idx) > 0) {
    # Columns to set to NA (excluding 'item', 'model', 'number', 'pvalue')
    cols_to_NA <- which(!names(out_items) %in% c("item", "model", "number", "pvalue"))
    out_items[bad_items_idx, cols_to_NA] <- NA
    if(con$verbose) message("\nWarning: Some items have 0 variance. Parameters set to NA.")
  }

  out_persons <- data.frame(
    id = 1:N,
    theta = round(final_theta, 3),
    se = round(final_se, 3),
    n_resp = rowSums(!is.na(raw_mat))
  )

  # Remove Row Names
  row.names(out_items) <- NULL
  row.names(out_persons) <- NULL

  n_par <- sum(ifelse(model %in% binary_models, ifelse(model=="Rasch",1, ifelse(model=="2PL",2,3)), 1 + (n_cats-1)))
  if(length(bad_items_idx)>0) n_par <- n_par - length(bad_items_idx)

  aic <- 2*n_par - 2*log_lik
  bic <- n_par*log(N) - 2*log_lik

  fit_df <- data.frame(
    Stat = c("LogLikelihood", "AIC", "BIC", "Iterations"),
    Value = c(round(log_lik, 2), round(aic, 2), round(bic, 2), iter)
  )

  # --- Final Verbose Message ---
  if(con$verbose) {
    cat("Finished All Estimation.\n")
    cat("------------------------------------------------\n")
  }

  return(list(
    item_params = out_items,
    person_params = out_persons,
    model_fit = fit_df,
    settings = list(method=method, models=table(model))
  ))
}
