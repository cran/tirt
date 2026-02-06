#' Joint Item Response Theory and Testlet Response Theory Estimation (Dichotomous & Polytomous)
#'
#' @description
#' Provides a unified marginal maximum likelihood estimation framework for
#' a broad class of item response theory and testlet response theory models. The function automatically detects data
#' structures to apply appropriate models, along with their testlet-effect extensions (Bradlow et al., 1999).
#'
#' @param data A \code{data.frame} or \code{matrix} containing item responses.
#'   Responses should be 0-indexed integers. Missing values should be coded as \code{NA}.
#' @param item_spec A \code{data.frame} providing item metadata. Must include
#'   columns \code{"item"} (matching \code{colnames(data)}) and \code{"model"}.
#'   Optionally includes a \code{"testlet"} column for TRT specifications.
#' @param method A character string specifying the estimation method. Currently
#'   supports \code{"EM"} (Expectation-Maximization). Defaults to \code{"EM"}.
#' @param control A \code{list} of control parameters for the estimation algorithm:
#'   \itemize{
#'     \item \code{max_iter}: Maximum number of EM iterations (default = 100).
#'     \item \code{converge_tol}: Convergence criterion for parameter change (default = 1e-4).
#'     \item \code{theta_range}: Numeric vector of length 2 specifying the integration
#'           grid bounds (default = c(-4, 4)).
#'     \item \code{quad_points}: Number of quadrature points (default = 21).
#'     \item \code{verbose}: Logical; if \code{TRUE}, prints progress to console.
#'     \item \code{fix_discrimination}: Logical; default=FALSE
#'   }
#'
#' @details
#' The estimation utilizes a robust Newton-Raphson update within the M-step. For
#' testlet models, dimension reduction is achieved through the integration of
#' the nuisance testlet effect (Li et al., 2006). The function automatically
#' corrects model specifications if the data levels (binary vs. polytomous)
#' do not align with the requested model string.
#'
#' @return A \code{list} containing three components:
#' \item{item_params}{A data frame of estimated item slopes (discrimination),
#'   difficulties/thresholds, and guessing parameters with associated standard errors.}
#' \item{person_params}{A data frame of EAP-based ability estimates (\eqn{\theta})
#'   and testlet effect estimates (\eqn{\gamma}).}
#' \item{model_fit}{A data frame containing Log-Likelihood, AIC, and BIC indices.}
#'
#' @references
#' Bradlow, E. T., Wainer, H., & Wang, X. (1999). A testlet response model for
#' multidimensionality in item response theory. \emph{Psychometrika, 64}(2), 147-168.
#'
#' Li, Y., Bolt, D. M., & Fu, J. (2006). A comparison of methods for estimating
#' secondary dimensions in testlet-based data. \emph{Applied Psychological
#' Measurement, 30}(3), 203-223.
#'
#' @importFrom stats dnorm qlogis var setNames
#' @examples
#'   # --- Example: Simulation (Binary + Poly + Testlets) ---
#'   set.seed(2025)
#'   N <- 100; J <- 20
#'
#'   # 1. Generate Parameters
#'   theta <- rnorm(N, 0, 1)
#'   gamma_1 <- rnorm(N, 0, 0.5) # Testlet 1 effect
#'   gamma_2 <- rnorm(N, 0, 0.6) # Testlet 2 effect
#'
#'   a_true <- runif(J, 0.8, 1.5)
#'   b_true <- seq(-1.5, 1.5, length.out = J)
#'
#'   resp_matrix <- matrix(NA, N, J)
#'   colnames(resp_matrix) <- paste0("Item_", 1:J)
#'
#'   # 2. Simulate Responses
#'   # Items 1-10: Binary Independent (Model: 2PL)
#'   for(j in 1:10) {
#'     p <- 1 / (1 + exp(-a_true[j] * (theta - b_true[j])))
#'     resp_matrix[,j] <- rbinom(N, 1, p)
#'   }
#'
#'   # Items 11-15: Poly Independent (Model: GRM)
#'   for(j in 11:15) {
#'     thresh <- sort(c(b_true[j] - 0.7, b_true[j] + 0.7))
#'     p1 <- 1 / (1 + exp(-a_true[j] * (theta - thresh[1])))
#'     p2 <- 1 / (1 + exp(-a_true[j] * (theta - thresh[2])))
#'     probs <- cbind(1-p1, p1-p2, p2)
#'     resp_matrix[,j] <- apply(probs, 1, function(p) sample(0:2, 1, prob=p))
#'   }
#'
#'   # Items 16-17: Binary Testlet 1 (Model: 2PLT)
#'   for(j in 16:17) {
#'     eff_theta <- theta + gamma_1
#'     p <- 1 / (1 + exp(-a_true[j] * (eff_theta - b_true[j])))
#'     resp_matrix[,j] <- rbinom(N, 1, p)
#'   }
#'
#'   # Items 18-20: Poly Testlet 2 (Model: GRT)
#'   for(j in 18:20) {
#'     eff_theta <- theta + gamma_2
#'     thresh <- sort(c(b_true[j] - 0.5, b_true[j] + 0.5))
#'     p1 <- 1 / (1 + exp(-a_true[j] * (eff_theta - thresh[1])))
#'     p2 <- 1 / (1 + exp(-a_true[j] * (eff_theta - thresh[2])))
#'     probs <- cbind(1-p1, p1-p2, p2)
#'     resp_matrix[,j] <- apply(probs, 1, function(p) sample(0:2, 1, prob=p))
#'   }
#'
#'   df_sim <- as.data.frame(resp_matrix)
#'
#'   # 3. Create Item Specification
#'   # STRICT naming: Independent=2PL/GRM, Testlet=2PLT/GRT
#'   spec <- data.frame(
#'     item = colnames(df_sim),
#'     model = c(rep("2PL", 10), rep("GRM", 5), rep("2PLT", 2), rep("GRT", 3)),
#'     testlet = c(rep(NA, 15), rep("T1", 2), rep("T2", 3)),
#'     stringsAsFactors = FALSE
#'   )
#'
#'   # 4. Run Estimation
#'   res <- irt_trt(df_sim, spec, method = "EM",
#'                  control = list(max_iter = 20, verbose = FALSE))
#'
#'   head(res$item_params)
#'   head(res$person_params)
#' @export
irt_trt <- function(data,
                          item_spec,
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
    verbose = TRUE,
    fix_discrimination = FALSE
  )
  con[names(control)] <- control

  # --- Input Validation ---
  if(!is.data.frame(data)) stop("Input 'data' must be a data frame.")
  raw_mat <- as.matrix(data)
  if(!is.numeric(raw_mat)) stop("Input data must be numeric.")

  N <- nrow(raw_mat)
  J <- ncol(raw_mat)
  item_names <- colnames(data)
  if(is.null(item_names)) item_names <- paste0("Item_", 1:J)
  colnames(raw_mat) <- item_names

  # Validate spec
  if(!is.data.frame(item_spec)) stop("'item_spec' must be a data.frame.")
  spec_indices <- match(item_names, item_spec$item)
  if(any(is.na(spec_indices))) stop("Error: Item names in data do not match 'item_spec'.")

  item_spec <- item_spec[spec_indices, ]
  mod_vec <- as.character(item_spec$model)

  # Testlet setup
  testlet_vec <- rep(NA, J)
  t_names_map <- NULL
  if("testlet" %in% names(item_spec)) {
    raw_t <- as.character(item_spec$testlet)
    raw_t[raw_t == "" | is.na(raw_t)] <- NA
    unique_t <- unique(raw_t[!is.na(raw_t)])
    if(length(unique_t) > 0) {
      t_names_map <- unique_t
      t_map <- setNames(1:length(unique_t), unique_t)
      testlet_vec <- as.integer(t_map[raw_t])
    }
  }
  n_testlets <- if(all(is.na(testlet_vec))) 0 else max(testlet_vec, na.rm=TRUE)

  # Shift Polytomous Data
  min_val <- min(raw_mat, na.rm=TRUE)
  if(is.finite(min_val) && min_val > 0) {
    if(con$verbose) cat("Note: Input responses shifted to start at 0 (internal requirement).\n")
    raw_mat <- raw_mat - min_val
  }
  n_cats <- apply(raw_mat, 2, function(x) max(x, na.rm=TRUE) + 1)

  valid_items_idx <- which(apply(raw_mat, 2, var, na.rm=TRUE) > 0)


  mismatched_items <- c()

  for(j in 1:J) {
    old_mod <- mod_vec[j]
    new_mod <- old_mod # default assumption
    cats <- n_cats[j]

    # -----------------------------------------------------
    # CASE 1: Data is Polytomous (cats > 2), but Model is Binary
    # -----------------------------------------------------
    if(cats > 2) {
      # Standard Binary -> Poly Equivalent
      if(old_mod %in% c("2PL", "3PL", "1PL")) new_mod <- "GPCM"
      else if(old_mod == "Rasch") new_mod <- "PCM"

      # Testlet Binary -> Poly Testlet Equivalent
      else if(old_mod %in% c("2PLT", "3PLT", "1PLT")) new_mod <- "GPCMT"
      else if(old_mod == "RaschT") new_mod <- "PCMT"
    }

    # -----------------------------------------------------
    # CASE 2: Data is Binary (cats == 2), but Model is Poly
    # -----------------------------------------------------
    else if(cats == 2) {
      # Poly -> Binary Equivalent
      if(old_mod %in% c("GPCM", "GRM")) new_mod <- "2PL"
      else if(old_mod == "PCM") new_mod <- "Rasch"

      # Poly Testlet -> Binary Testlet Equivalent
      else if(old_mod %in% c("GPCMT", "GRT")) new_mod <- "2PLT"
      else if(old_mod == "PCMT") new_mod <- "RaschT"
    }

    # Apply Correction
    if(new_mod != old_mod) {
      mod_vec[j] <- new_mod
      mismatched_items <- c(mismatched_items, sprintf("%s (Data: %d cats) | Spec: %s -> Corrected: %s",
                                                      item_names[j], cats, old_mod, new_mod))
    }
  }


  # --- User-Friendly Start Message ---
  if(con$verbose) {
    cat("\n================================================\n")
    cat("        UNIVERSAL IRT/TRT ESTIMATION            \n")
    cat("================================================\n")

    # Summary of configurations
    n_binary <- sum(n_cats == 2)
    n_poly   <- sum(n_cats > 2)
    n_indep  <- sum(is.na(testlet_vec))
    n_testlet_items <- sum(!is.na(testlet_vec))

    cat("Data Summary:\n")
    cat(sprintf("  - Total Items: %d\n", J))
    cat(sprintf("  - Binary Items: %d\n", n_binary))
    cat(sprintf("  - Polytomous Items: %d\n", n_poly))
    cat("\nStructure:\n")
    cat(sprintf("  - Independent Items: %d\n", n_indep))
    cat(sprintf("  - Testlet Items: %d (grouped into %d testlets)\n", n_testlet_items, n_testlets))

    if(n_testlets > 0) {
      cat("\nNOTE on Testlets:\n")
      cat("  Any item with a value in the 'testlet' string within the item_spec statement will be estimated\n")
      cat("  as a Testlet item (Theta + Gamma), regardless of the model string.\n")
      cat("  (e.g., model='GPCM' + testlet='T1' -> Testlet GPCT).\n")
    }

    # --- Print Warnings for Mismatches ---
    if(length(mismatched_items) > 0) {
      cat("\n[!] WARNING: MODEL SPECIFICATION MISMATCH\n")
      cat("    The following items had models inconsistent with their data structure.\n")
      cat("    They have been automatically corrected to the appropriate model:\n")
      cat(paste("    -", mismatched_items, collapse="\n"))
      cat("\n")
    }

    cat("------------------------------------------------\n")
  }

  # --- Initialize Parameters ---
  a <- rep(1.0, J)
  b_list <- vector("list", J)
  g <- rep(0.0, J)

  # Initial Estimates (CTT)
  p_vals <- colMeans(raw_mat, na.rm=TRUE) / (n_cats - 1)
  p_vals <- pmin(pmax(p_vals, 0.01), 0.99)

  for(j in 1:J) {
    mod <- mod_vec[j]
    K <- n_cats[j]
    if(mod %in% c("Rasch", "RaschT", "PCM", "PCMT")) a[j] <- 1.0 else a[j] <- 1.0

    y <- raw_mat[,j]; y <- y[!is.na(y)]
    if(length(y) < 2) { b_list[[j]] <- rep(0, max(1, K-1)); next }

    if(K == 2) {
      b_list[[j]] <- -log(p_vals[j] / (1 - p_vals[j]))
      if(mod %in% c("3PL", "3PLT")) g[j] <- 0.2
    } else {
      props <- table(factor(y, levels=0:(K-1))) / length(y)
      cum <- cumsum(props)[1:(K-1)]; cum <- pmin(pmax(cum, 0.01), 0.99)
      thresh <- -qlogis(cum)
      if(mod %in% c("GRM", "GRT")) b_list[[j]] <- sort(thresh)
      else b_list[[j]] <- seq(-0.5, 0.5, length.out=K-1)
    }
  }

  sigma_gamma <- rep(0.5, max(1, n_testlets))

  # MLE storage (used if method=MLE)
  theta_mle <- as.vector(scale(rowSums(raw_mat, na.rm=TRUE)))
  theta_mle[is.na(theta_mle)] <- 0
  gamma_mle <- matrix(0, N, max(1, n_testlets))

  # --- Helper: Universal Probability ---
  plogis_c <- function(x) 1 / (1 + exp(-x))

  get_P_univ <- function(th_eff, a_val, b_vec, g_val, mod, n_cat) {
    n_q <- length(th_eff)
    if(n_cat == 2) {
      z <- a_val * (th_eff - b_vec[1])
      z <- pmin(pmax(z, -30), 30)
      p <- g_val + (1 - g_val) * plogis_c(z)
      return(cbind(1-p, p))
    }
    if(mod %in% c("GRM", "GRT")) {
      probs <- matrix(0, n_q, n_cat)
      P_s <- matrix(0, n_q, n_cat+1); P_s[,1] <- 1
      for(k in 1:(n_cat-1)) {
        z <- a_val * (th_eff - b_vec[k]); z <- pmin(pmax(z, -30), 30)
        P_s[, k+1] <- plogis_c(z)
      }
      for(k in 1:n_cat) probs[,k] <- P_s[,k] - P_s[,k+1]
      return(pmax(probs, 1e-15))
    }
    # GPCM/PCM
    numer <- matrix(0, n_q, n_cat)
    curr <- 0
    for(k in 1:(n_cat-1)) {
      curr <- curr + a_val * (th_eff - b_vec[k])
      numer[, k+1] <- curr
    }
    mx <- apply(numer, 1, max)
    exps <- exp(numer - mx)
    return(exps / rowSums(exps))
  }

  # --- Helper: Robust Newton-Raphson Item Update ---
  solve_item <- function(counts, nodes, ca, cb, cg, mod, n_cat, fix_d) {
    start_p <- c()
    if(!fix_d) start_p <- c(start_p, ca)
    start_p <- c(start_p, cb)
    if(mod %in% c("3PL", "3PLT")) start_p <- c(start_p, cg)

    n_par <- length(start_p)

    # Initialize Hessian to avoid scoping crash
    he <- matrix(0, n_par, n_par)
    diag(he) <- -1

    ll_fn <- function(p) {
      idx <- 1
      la <- if(fix_d) 1.0 else { v<-p[idx]; idx<-idx+1; v }
      lb <- p[idx:(idx+length(cb)-1)]; idx<-idx+length(cb)
      lg <- if(mod %in% c("3PL", "3PLT")) p[idx] else 0

      if(la < 0.01 || la > 6.0) return(-1e20)
      if(mod %in% c("3PL", "3PLT") && (lg < 0 || lg > 0.45)) return(-1e20)
      if(mod %in% c("GRM", "GRT") && is.unsorted(lb)) return(-1e20)

      pr <- get_P_univ(nodes, la, lb, lg, mod, n_cat)
      # pmax to prevent log(0)
      sum(counts * log(pmax(pr, 1e-100)))
    }

    curr <- start_p
    opt_fail <- FALSE

    for(iter in 1:con$nr_max_iter) {
      f0 <- ll_fn(curr)
      if(f0 <= -1e19) { opt_fail <- TRUE; break }

      gr <- numeric(n_par)
      he <- matrix(0, n_par, n_par)
      h <- 1e-4

      for(i in 1:n_par) {
        p1 <- curr; p1[i] <- p1[i] + h
        gr[i] <- (ll_fn(p1) - f0)/h
      }
      for(i in 1:n_par) {
        p1 <- curr; p1[i] <- p1[i] + h; p2 <- curr; p2[i] <- p2[i] - h
        he[i,i] <- (ll_fn(p1) - 2*f0 + ll_fn(p2))/h^2
      }
      diag(he) <- diag(he) - 1e-3

      step <- try(solve(he, gr), silent=TRUE)
      if(inherits(step, "try-error")) step <- gr * 0.01
      step <- pmax(pmin(step, 1.0), -1.0)

      curr <- curr - step * con$nr_damp
      if(max(abs(step)) < 1e-4) break
    }

    se <- tryCatch({ sqrt(diag(solve(-he + diag(1e-5, n_par)))) }, error=function(e) rep(NA, n_par))

    if(opt_fail) {
      return(list(se=rep(NA, n_par), a=ca, b=cb, g=cg))
    }

    res <- list(se=se); idx <- 1
    if(!fix_d) { res$a <- curr[idx]; idx<-idx+1 } else res$a <- 1.0
    res$b <- curr[idx:(idx+length(cb)-1)]; idx<-idx+length(cb)
    if(mod %in% c("3PL","3PLT")) res$g <- curr[idx] else res$g <- 0
    return(res)
  }

  # --- Main Estimation Loop ---
  if(con$verbose) {
    cat(sprintf("Data: %d Persons, %d Items, %d Testlets\n", N, J, n_testlets))
    cat("------------------------------------------------\n")
    cat("Starting Estimation Loop...\n")
  }

  nodes <- seq(con$theta_range[1], con$theta_range[2], length.out=con$quad_points)
  wts <- dnorm(nodes); wts <- wts/sum(wts)

  log_lik <- 0
  se_list <- vector("list", J)
  Post <- matrix(1/length(nodes), N, length(nodes)) # Init posterior

  is_converged <- FALSE

  for(iter in 1:con$max_iter) {

    # E-STEP
    L_tot <- matrix(1, N, length(nodes))

    # Independent Items
    for(j in valid_items_idx) {
      if(!is.na(testlet_vec[j])) next
      P <- get_P_univ(nodes, a[j], b_list[[j]], g[j], mod_vec[j], n_cats[j])
      y <- raw_mat[,j]; valid <- !is.na(y)
      if(any(valid)) {
        for(q in 1:length(nodes)) L_tot[valid, q] <- L_tot[valid, q] * P[q, y[valid]+1]
      }
    }

    # Testlet Items (Dim Reduction for EM)
    if(n_testlets > 0) {
      g_grid <- seq(-3, 3, length.out=9)
      gw <- dnorm(g_grid); gw <- gw/sum(gw)

      for(t in 1:n_testlets) {
        t_idx <- which(testlet_vec == t)
        L_t <- matrix(0, N, length(nodes))

        for(q in 1:length(nodes)) {
          eff_t <- nodes[q] + sqrt(sigma_gamma[t]) * g_grid
          L_sub <- matrix(1, N, length(g_grid))
          for(j in t_idx) {
            P <- get_P_univ(eff_t, a[j], b_list[[j]], g[j], mod_vec[j], n_cats[j])
            y <- raw_mat[,j]; valid <- !is.na(y)
            if(any(valid)) {
              for(k in 1:length(g_grid)) L_sub[valid, k] <- L_sub[valid, k] * P[k, y[valid]+1]
            }
          }
          L_t[, q] <- as.vector(L_sub %*% gw)
        }
        L_tot <- L_tot * L_t
      }
    }

    Ly <- as.vector(L_tot %*% wts)
    log_lik <- sum(log(pmax(Ly, 1e-100)))
    Post <- (L_tot * matrix(wts, N, length(nodes), byrow=TRUE)) / Ly

    # M-STEP
    max_diff <- 0
    for(j in valid_items_idx) {
      K <- n_cats[j]
      counts <- matrix(0, length(nodes), K)
      y <- raw_mat[,j]
      valid <- !is.na(y)
      for(k in 1:K) {
        idx <- which(valid & y == (k-1))
        if(length(idx)>0) counts[,k] <- colSums(Post[idx, , drop=FALSE])
      }

      old <- c(a[j], b_list[[j]], g[j])
      res <- solve_item(counts, nodes, a[j], b_list[[j]], g[j], mod_vec[j], K, con$fix_discrimination)
      a[j] <- res$a; b_list[[j]] <- res$b; g[j] <- res$g; se_list[[j]] <- res$se

      max_diff <- max(max_diff, abs(c(a[j], b_list[[j]], g[j]) - old))
    }

    if(con$verbose) cat(sprintf("\rIter %3d | LogLik: %10.2f | Max Change: %.5f    ", iter, log_lik, max_diff))

    if(max_diff < con$converge_tol) {
      if(con$verbose) cat("\n\n>>> Convergence Confirmation: Model Converged!\n")
      is_converged <- TRUE
      break
    }
  }

  if(!is_converged && con$verbose) {
    cat(sprintf("\n\n>>> NOTE: Stopped because max Iteration Reached at %d.\n", con$max_iter))
    cat("    You can increase 'max_iter' in control=list(max_iter=...) if needed.\n")
  }

  # --- Post-Hoc: Person Parameter Estimation ---
  if(con$verbose) {
    cat("------------------------------------------------\n")
    cat("Estimating Final Person Parameters (Theta)...\n")
    cat(sprintf("Using Theta range: [%.1f, %.1f]. Adjust via control=list(theta_range=c(...)).\n", con$theta_range[1], con$theta_range[2]))
  }

  theta_final <- as.vector(Post %*% nodes)
  se_final <- sqrt(abs(as.vector(Post %*% nodes^2) - theta_final^2))

  # --- Post-Hoc: Testlet Effect Estimation (Gamma) ---
  gamma_final <- matrix(NA, N, max(1, n_testlets))
  gamma_se_final <- matrix(NA, N, max(1, n_testlets))

  if(n_testlets > 0) {
    if(con$verbose) cat("Estimating Testlet Effects (Gamma) per person...\n")

    g_grid <- seq(-3, 3, length.out=15) # Finer grid for EAP
    gw <- dnorm(g_grid); gw <- gw/sum(gw)

    for(t in 1:n_testlets) {
      # For EAP(Gamma): We integrate Gamma posterior.
      # P(gamma | y) \propto P(y_t | theta, gamma) * P(theta | y_others) * P(gamma)
      # Approx: use P(theta | all_y) as weight, though slightly circular, is standard approx in unidi code.
      # Better: Joint EAP over Theta/Gamma grid weighted by Posterior Theta.

      t_idx <- which(testlet_vec == t)

      # For each person, we construct a marginal distribution over Gamma
      # Joint(q, k) = Post_Theta[i, q] * P(y_t | theta_q, gamma_k) * P(gamma_k)
      # Note: We divide by L_t(theta_q) to remove double counting if strictly needed,
      # but standard scoring often just uses the joint expectation.

      for(i in 1:N) {
        valid_resp <- !is.na(raw_mat[i, t_idx])
        if(!any(valid_resp)) { gamma_final[i, t] <- 0; gamma_se_final[i, t] <- NA; next }

        # Calculate Likelihood of testlet items on the 2D grid (Theta x Gamma)
        # We perform summation over Theta index (q) weighted by Post[i, q]

        gam_post_num <- numeric(length(g_grid))

        # Pre-calc item probs for this person's responses
        # To speed up: assume Theta is fixed at EAP Theta? (Faster)
        # Or Integrate? (Better). Let's Integrate.

        for(k in 1:length(g_grid)) {
          g_val <- sqrt(sigma_gamma[t]) * g_grid[k]

          # Likelihood of this gamma given theta_q
          L_g_given_th <- numeric(length(nodes))

          # Vectorized P calc
          # Eff Theta Matrix [Nodes x 1]
          eff_mat <- nodes + g_val

          # Accumulate log-lik for items
          log_l <- rep(0, length(nodes))
          for(j_sub in t_idx) {
            if(is.na(raw_mat[i, j_sub])) next
            p_vec <- get_P_univ(eff_mat, a[j_sub], b_list[[j_sub]], g[j_sub], mod_vec[j_sub], n_cats[j_sub])
            log_l <- log_l + log(p_vec[, raw_mat[i, j_sub]+1])
          }
          lik_vec <- exp(log_l)

          # Weighted sum over Theta Posterior
          # P(gamma_k) ~ Sum_q ( P(y_t | th_q, gam_k) * Post(th_q) ) * Prior(gam_k)
          gam_post_num[k] <- sum(lik_vec * Post[i, ]) * gw[k]
        }

        denom <- sum(gam_post_num)
        if(denom == 0) denom <- 1e-10
        gam_dist <- gam_post_num / denom

        # EAP
        real_g_vals <- sqrt(sigma_gamma[t]) * g_grid
        gamma_final[i, t] <- sum(real_g_vals * gam_dist)
        gamma_se_final[i, t] <- sqrt(sum((real_g_vals - gamma_final[i, t])^2 * gam_dist))
      }
    }
  }

  if(con$verbose) cat("All parameters estimated. Formatting results...\n")

  # --- Output Construction ---

  # 1. Item Params
  get_se <- function(j, type, k=1) {
    v <- se_list[[j]]
    if(is.null(v) || any(is.na(v))) return(NA)
    has_a <- !con$fix_discrimination
    has_g <- mod_vec[j] %in% c("3PL", "3PLT")
    n_b <- length(b_list[[j]])
    idx_a <- if(has_a) 1 else 0
    idx_b_start <- idx_a + 1
    idx_g <- if(has_g) idx_b_start + n_b else 0

    if(type == "a") return(if(has_a) v[1] else NA)
    if(type == "g") return(if(has_g) v[idx_g] else NA)
    if(type == "b") {
      if(k > n_b) return(NA)
      return(v[idx_b_start + k - 1])
    }
    return(NA)
  }

  df_items <- data.frame(
    item = item_names,
    model = mod_vec,
    testlet = ifelse(is.na(testlet_vec), NA,
                     if(!is.null(t_names_map)) t_names_map[testlet_vec] else paste0("T", testlet_vec)),
    discrimination = round(a, 3),
    discrimination_se = round(sapply(1:J, get_se, type="a"), 3),

    # --- Binary Items use Difficulty, Poly uses NA ---
    difficulty = round(sapply(1:J, function(j) if(n_cats[j] == 2) b_list[[j]][1] else NA), 3),
    difficulty_se = round(sapply(1:J, function(j) if(n_cats[j] == 2) get_se(j, type="b", k=1) else NA), 3),
    # -----------------------------------------------------------------

    guessing = round(g, 3),
    guessing_se = round(sapply(1:J, get_se, type="g"), 3),
    stringsAsFactors = FALSE
  )

  # Add thresholds (ONLY for Polytomous, to keep clean separation)
  max_thresh <- max(sapply(b_list, length))
  for(k in 1:max_thresh) {
    # Extract value only if item is Polytomous
    val <- sapply(1:J, function(j) {
      if(n_cats[j] > 2 && length(b_list[[j]]) >= k) b_list[[j]][k] else NA
    })
    se_val <- sapply(1:J, function(j) {
      if(n_cats[j] > 2 && length(b_list[[j]]) >= k) get_se(j, type="b", k=k) else NA
    })

    df_items[[paste0("step/threshold_", k)]] <- round(val, 3)
    df_items[[paste0("step/threshold_", k, "_se")]] <- round(se_val, 3)
  }

  # Add number and pvalue at the end
  df_items$number <- colSums(!is.na(raw_mat))
  df_items$pvalue <- round(colMeans(raw_mat, na.rm=TRUE)/(n_cats-1), 3)

  # Remove row names
  row.names(df_items) <- NULL

  # 2. Person Params
  df_persons <- data.frame(
    n_response = rowSums(!is.na(raw_mat)),
    ability = round(theta_final, 3),
    ability_se = round(se_final, 3)
  )

  if(n_testlets > 0) {
    for(t in 1:n_testlets) {
      t_name <- if(!is.null(t_names_map)) t_names_map[t] else paste0("testlet_", t)
      df_persons[[t_name]] <- round(gamma_final[, t], 3)
      df_persons[[paste0(t_name, "_se")]] <- round(gamma_se_final[, t], 3)
    }
  }
  row.names(df_persons) <- NULL

  # 3. Model Fit
  n_par <- length(valid_items_idx) * 2 # Approximation
  df_fit <- data.frame(
    Index=c("LogLik","AIC","BIC", "Iterations"),
    Value=c(round(log_lik,3), round(2*n_par - 2*log_lik,3), round(n_par*log(N) - 2*log_lik,3), iter)
  )

  if(con$verbose) {
    cat("All completed. You can view the results now.\n")
    cat("================================================\n")
  }

  return(list(item_params=df_items, person_params=df_persons, model_fit=df_fit))
}
