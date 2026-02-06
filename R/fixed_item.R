#' Fixed Item Calibration
#'
#' @description
#' Estimates unknown item parameters using Marginal Maximum Likelihood via Expectation-Maximization Algorithm.
#' Uses a custom Bounded Newton-Raphson solver.
#' Supports mixed-format data containing dichotomous and polytomous responses
#'
#' @param response_df A data.frame of responses. Rows=Students, Cols=Items.
#'   Data MUST be from 0-indexed (0, 1, 2...).
#' @param item_params_df A data.frame of known parameters.
#'   Required: "item", "model".
#' @param control A \code{list} of control parameters for the estimation algorithm:
#'   \itemize{
#'     \item \code{max_iter}: Maximum number of EM iterations (default = 50).
#'     \item \code{conv_crit}: Convergence criterion for parameter change (default = 0.005).
#'     \item \code{verbose}: Logical; if \code{TRUE}, prints progress to console.
#'   }
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{item_params}: Estimated parameters for unknown items.
#'   \item \code{person_params}: Estimated person parameters.
#'   \item \code{model_fit}: A data frame containing number off estimated parameters and fit statistics such as Akaike’s Information Criterion (AIC), the Bayesian Information Criterion (BIC), and Log-Likelihood.
#' }
#'
#' @examples
#' # 1. TOY EXAMPLE
#' # ===========================================================================
#' set.seed(123)
#' # Create a very small dataset (N=50, J=4)
#' N_toy <- 50
#' df_toy <- data.frame(
#'   I1 = rbinom(N_toy, 1, 0.5), I2 = rbinom(N_toy, 1, 0.6), # Known items
#'   U1 = rbinom(N_toy, 1, 0.5), U2 = rbinom(N_toy, 1, 0.4)  # Unknown items
#' )
#'
#' # Define the "Known" parameters for I1 and I2
#' known_params <- data.frame(
#'   item = c("I1", "I2"),
#'   model = c("2PL", "2PL"),
#'   a = c(1.0, 1.2),
#'   b = c(-0.5, 0.5)
#' )
#'
#' # Run Fixed Item Calibration with very low iterations
#' fit_toy <- fixed_item(df_toy, known_params, control=list(max_iter=2, verbose=FALSE))
#' print(head(fit_toy$item_params))
#'   \donttest{
#'   # --- Example 2: Simulation ---
#'   set.seed(123)
#'   N <- 500
#'   true_theta <- rnorm(N, 0, 1)
#'
#'   # 1. Simulation Helpers
#'   sim_2pl <- function(theta, a, b) {
#'     p <- 1 / (1 + exp(-1.7 * a * (theta - b)))
#'     rbinom(N, 1, p)
#'   }
#'   sim_poly <- function(theta, a, steps) {
#'     n_cat <- length(steps) + 1
#'     probs <- matrix(0, length(theta), n_cat)
#'     for(k in 1:n_cat) {
#'       score <- k - 1
#'       if(score == 0) num <- 0
#'       else num <- a * (score * theta - sum(steps[1:score]))
#'       probs[, k] <- exp(num)
#'     }
#'     probs <- probs / rowSums(probs)
#'     apply(probs, 1, function(x) sample(0:(n_cat-1), 1, prob=x))
#'   }
#'
#'   # 2. Generate Data (Mixed Known/Unknown Items)
#'   # Items 1-5: Known Binary (2PL)
#'   # Items 6-10: Unknown Binary (2PL)
#'   # Items 11-12: Known Poly (GPCM)
#'   # Items 13-15: Unknown Poly (GPCM)
#'
#'   resp_mat <- matrix(NA, N, 15)
#'   colnames(resp_mat) <- paste0("Item_", 1:15)
#'
#'   # Known Binary Parameters
#'   a_bin <- c(1.0, 1.2, 0.9, 1.1, 0.8)
#'   b_bin <- c(-1, -0.5, 0, 0.5, 1)
#'
#'   for(i in 1:5) resp_mat[,i] <- sim_2pl(true_theta, a_bin[i], b_bin[i])
#'   for(i in 6:10) resp_mat[,i] <- sim_2pl(true_theta, runif(1,0.8,1.2), rnorm(1))
#'
#'   # Known Poly Parameters
#'   a_poly <- c(1.0, 0.8)
#'   d_poly <- list(c(-1, 1), c(-0.5, 0.5))
#'
#'   resp_mat[,11] <- sim_poly(true_theta, a_poly[1], d_poly[[1]])
#'   resp_mat[,12] <- sim_poly(true_theta, a_poly[2], d_poly[[2]])
#'   for(i in 13:15) resp_mat[,i] <- sim_poly(true_theta, 1.0, c(-0.5, 0.5))
#'
#'   df_resp <- as.data.frame(resp_mat)
#'
#'   # 3. Create 'Known Parameters' Dataframe
#'   # This tells the function: "Fix these, Estimate the rest"
#'   known_df <- data.frame(
#'     item = c(paste0("Item_", 1:5), "Item_11", "Item_12"),
#'     model = c(rep("2PL", 5), rep("GPCM", 2)),
#'     a = c(a_bin, a_poly),
#'     b = c(b_bin, NA, NA),     # Binary difficulty
#'     step_1 = c(rep(NA, 5), -1, -0.5), # Poly steps
#'     step_2 = c(rep(NA, 5),  1,  0.5),
#'     stringsAsFactors = FALSE
#'   )
#'
#'   # 4. Run Estimation
#'   res <- fixed_item(df_resp, known_df, control=list(max_iter=20))
#'
#'   # View Results
#'   # Notice Items 1-5 and 11-12 have Status "Fixed"
#'   head(res$item_params, 12)
#'
#'   # --- Example 2: With Package Data ---
#'   data("ela1", package = "tirt")
#'
#'   # Let's treat the first 5 items as "Known" with arbitrary parameters
#'   # just to demonstrate syntax.
#'   df_real <- ela1[, 1:20]
#'
#'   known_real <- data.frame(
#'     item = paste0("Q", 1:5),
#'     model = "2PL",
#'     a = 1.0,
#'     b = seq(-1, 1, length.out=5)
#'   )
#'
#'   # Ideally, column names in df_real should match 'item' column in known_real
#'   colnames(df_real)[1:5] <- paste0("Q", 1:5)
#'
#'   real_res <- fixed_item(df_real, known_real, control=list(max_iter=10))
#'   head(real_res$item_params)
#'   }
#' @export
fixed_item <- function(response_df, item_params_df, control = list()) {

  # --- 0. Control Parameters ---
  con <- list(max_iter = 50, conv_crit = 0.005, verbose = TRUE)
  con[names(control)] <- control

  # --- 1. INTERNAL MATH FUNCTIONS (The Engine) ---

  # 1.1 Probability Engine
  get_prob_matrix <- function(theta, a, b_vec, c_param, model) {
    n_q <- length(theta)

    # Binary Models
    if (model %in% c("Rasch", "2PL", "3PL")) {
      D <- 1.7
      lin <- D * a * (theta - b_vec[1])
      # Logistic function
      p1 <- 1 / (1 + exp(-lin))
      # 3PL Guessing
      if (model == "3PL") p1 <- c_param + (1 - c_param) * p1
      p0 <- 1 - p1
      return(cbind(p0, p1))
    }

    # GPCM / PCM (Divide by total)
    if (model %in% c("GPCM", "PCM")) {
      n_cat <- length(b_vec) + 1
      numerators <- matrix(0, nrow = n_q, ncol = n_cat)
      current_sum <- 0
      numerators[, 1] <- 0
      for (k in 2:n_cat) {
        term <- a * (theta - b_vec[k-1])
        current_sum <- current_sum + term
        numerators[, k] <- current_sum
      }
      exps <- exp(numerators)
      # Row sums safety to prevent overflow/NaN
      rs <- rowSums(exps)
      rs[rs == 0] <- 1
      probs <- exps / rs
      return(probs)
    }

    # GRM (Difference of cumulative)
    if (model == "GRM") {
      n_cat <- length(b_vec) + 1
      D <- 1.7
      p_star <- matrix(0, nrow = n_q, ncol = n_cat + 1)
      p_star[, 1] <- 1.0; p_star[, n_cat + 1] <- 0.0

      for (k in 1:(n_cat - 1)) {
        lin <- D * a * (theta - b_vec[k])
        p_star[, k+1] <- 1 / (1 + exp(-lin))
      }
      probs <- matrix(0, nrow = n_q, ncol = n_cat)
      for (k in 1:n_cat) probs[, k] <- p_star[, k] - p_star[, k+1]
      probs <- pmax(probs, 1e-10) # Floor small probs
      probs <- probs / rowSums(probs)
      return(probs)
    }
    stop(paste("Unknown Model:", model))
  }

  # 1.2 Custom Bounded Newton-Raphson Solver
  # Used for both M-Step (Item Params) and Person Scoring
  custom_solver <- function(start_params, fn_minimize, lower_b, upper_b) {
    p <- start_params
    n_p <- length(p)

    # Settings for the inner loop
    max_inner <- 10
    tol <- 1e-4
    h <- 1e-4 # Step size for derivatives

    final_hessian <- matrix(NA, n_p, n_p)

    for (iter in 1:max_inner) {
      val <- fn_minimize(p)

      grad <- numeric(n_p)
      hess <- numeric(n_p) # Diagonal approx for speed & stability

      # Compute Numerical Gradient & Diagonal Hessian
      for (k in 1:n_p) {
        p_plus <- p; p_plus[k] <- p_plus[k] + h
        p_minus <- p; p_minus[k] <- p_minus[k] - h

        v_plus <- fn_minimize(p_plus)
        v_minus <- fn_minimize(p_minus)

        # Central Difference Gradient
        grad[k] <- (v_plus - v_minus) / (2 * h)

        # Second Derivative (Hessian diag)
        h_val <- (v_plus - 2 * val + v_minus) / (h^2)
        hess[k] <- h_val
      }

      # Update Rule: p_new = p - LearningRate * (Grad / Hessian)
      # We ensure Hessian is positive (curvature checks)
      hess <- pmax(hess, 1e-6)

      update <- grad / hess

      # Damping: Don't jump too far
      update <- pmax(pmin(update, 1.0), -1.0)

      p_new <- p - update

      # Apply Bounds (Clamping)
      p_new <- pmax(pmin(p_new, upper_b), lower_b)

      # Check convergence
      if (max(abs(p - p_new)) < tol) {
        p <- p_new
        final_hessian <- diag(hess, nrow=n_p) # Store diagonal matrix
        break
      }
      p <- p_new

      # Final Hessian update on last loop
      if (iter == max_inner) final_hessian <- diag(hess, nrow=n_p)
    }

    return(list(par = p, hessian = final_hessian))
  }

  # 1.3 Invert Hessian for SE
  get_se <- function(hess_mat) {
    tryCatch({
      # Solve inversion with ridge for stability
      inv <- solve(hess_mat + diag(1e-8, nrow(hess_mat)))
      se <- sqrt(diag(inv))
      se[is.nan(se)] <- NA
      return(se)
    }, error = function(e) return(rep(NA, nrow(hess_mat))))
  }

  # --- 2. Input Setup ---

  if(con$verbose) {
    message("---------------------------------------------------------")
    message("Starting Fixed Item Calibration...")
    message("---------------------------------------------------------")
  }

  items <- colnames(response_df)
  resp_mat <- as.matrix(response_df)

  min_val <- min(resp_mat, na.rm = TRUE)
  if (min_val < 0) stop("Data must be 0-indexed (non-negative).")
  if (min_val > 0) stop("Data must start at 0 (e.g., 0/1).")

  known_items <- intersect(items, item_params_df$item)
  unknown_items <- setdiff(items, known_items)

  if (length(known_items) == 0) stop("No known items found.")

  # --- Intelligent Model Detection ---
  known_models <- unique(item_params_df$model)

  # Binary Hierarchy: 3PL > 2PL > Rasch
  if ("3PL" %in% known_models) {
    default_binary <- "3PL"
  } else if ("2PL" %in% known_models) {
    default_binary <- "2PL"
  } else {
    default_binary <- "Rasch"
  }

  # Poly Hierarchy: GRM > PCM > GPCM
  if ("GRM" %in% known_models) {
    default_poly <- "GRM"
  } else if ("PCM" %in% known_models) {
    default_poly <- "PCM"
  } else {
    default_poly <- "GPCM"
  }

  if(con$verbose) {
    message(paste("Inferred Unknown Binary:", default_binary))
    message(paste("Inferred Unknown Poly:  ", default_poly))
  }

  # --- 3. Initialize Parameters ---

  item_list <- list()
  total_est_params <- 0

  # 3.1 Load Known
  for (itm in known_items) {
    row <- item_params_df[item_params_df$item == itm, ]
    mod <- row$model[1]

    a <- if(!is.null(row$a) && !is.na(row$a)) row$a else 1.0
    c_p <- if(!is.null(row$c) && !is.na(row$c)) row$c else 0.0

    b_cols <- grep("^(b|d\\d+|step_\\d+)$", names(row), value=TRUE)
    b_cols <- b_cols[order(nchar(b_cols), b_cols)]
    b_vec <- as.numeric(row[b_cols])
    b_vec <- b_vec[!is.na(b_vec)]
    if (length(b_vec) == 0 && "b" %in% names(row)) b_vec <- row$b

    item_list[[itm]] <- list(
      type = "fixed", model = mod,
      a = a, b = b_vec, c = c_p,
      n_cat = length(b_vec) + 1,
      se_a=NA, se_b=rep(NA, length(b_vec)), se_c=NA
    )
    if(mod %in% c("Rasch","2PL","3PL")) item_list[[itm]]$n_cat <- 2
  }

  # 3.2 Initialize Unknown
  for (itm in unknown_items) {
    obs <- na.omit(resp_mat[, itm])
    max_k <- max(obs)
    n_cat <- max_k + 1
    mod_use <- if (max_k == 1) default_binary else default_poly

    # Initial guesses
    prop <- mean(obs)/max_k
    prop <- min(max(prop,0.01),0.99)
    logit_p <- -log(prop/(1-prop))

    a_init <- 1.0
    c_init <- 0.0
    b_init <- logit_p

    if (mod_use == "3PL") c_init <- 0.1
    if (!mod_use %in% c("Rasch","2PL","3PL")) {
      b_init <- seq(logit_p-0.5, logit_p+0.5, length.out=max_k)
    }

    # Count parameters
    n_p <- length(b_init)
    if (mod_use %in% c("2PL","3PL","GPCM","GRM")) n_p <- n_p + 1
    if (mod_use == "3PL") n_p <- n_p + 1
    total_est_params <- total_est_params + n_p

    item_list[[itm]] <- list(
      type = "estimated", model = mod_use,
      a = a_init, b = b_init, c = c_init,
      n_cat = n_cat,
      se_a=NA, se_b=rep(NA, length(b_init)), se_c=NA
    )
  }

  # --- 4. EM Algorithm ---

  nodes <- seq(-6, 6, length.out = 41)
  weights <- dnorm(nodes); weights <- weights / sum(weights)

  cycle <- 0
  converged <- FALSE
  last_ll <- -Inf

  while(!converged && cycle < con$max_iter) {
    cycle <- cycle + 1

    # --- E-Step ---
    posterior <- matrix(rep(weights, nrow(resp_mat)), nrow = nrow(resp_mat), byrow = TRUE)

    for (itm in items) {
      obj <- item_list[[itm]]
      probs <- get_prob_matrix(nodes, obj$a, obj$b, obj$c, obj$model)

      obs <- resp_mat[, itm]
      valid <- which(!is.na(obs))
      if (length(valid) == 0) next

      idx <- obs[valid] + 1
      lik_mat <- t(probs[, idx])
      posterior[valid, ] <- posterior[valid, ] * lik_mat
    }

    row_sums <- rowSums(posterior)
    row_sums[row_sums == 0] <- 1e-10
    posterior <- posterior / row_sums
    curr_ll <- sum(log(row_sums))

    # --- M-Step ---
    for (itm in unknown_items) {
      obj <- item_list[[itm]]
      obs <- resp_mat[, itm]
      valid <- which(!is.na(obs))
      r_qk <- matrix(0, nrow = length(nodes), ncol = obj$n_cat)

      if (length(valid) > 0) {
        clean_resp <- obs[valid]
        for (k in 0:(obj$n_cat - 1)) {
          is_k <- which(clean_resp == k)
          if (length(is_k) > 0) {
            subset_post <- posterior[valid[is_k], , drop=FALSE]
            r_qk[, k+1] <- colSums(subset_post)
          }
        }
      }

      # 1. Define Function to Minimize (Negative Log Likelihood)
      fn_neg_ll <- function(p) {
        t_a <- obj$a; t_b <- obj$b; t_c <- obj$c
        idx <- 1

        if (obj$model %in% c("2PL","3PL","GPCM","GRM")) { t_a <- p[idx]; idx <- idx+1 }
        if (obj$model == "3PL") { t_c <- p[idx]; idx <- idx+1 }
        n_b <- length(obj$b)
        t_b <- p[idx:(idx+n_b-1)]

        pr <- get_prob_matrix(nodes, t_a, t_b, t_c, obj$model)
        pr <- pmax(pr, 1e-10)
        return(-sum(r_qk * log(pr)))
      }

      # 2. Build Params and Bounds
      par_init <- c()
      lower_b <- c()
      upper_b <- c()

      if (obj$model %in% c("2PL","3PL","GPCM","GRM")) {
        par_init <- c(par_init, obj$a)
        lower_b <- c(lower_b, 0.01); upper_b <- c(upper_b, 8.0)
      }
      if (obj$model == "3PL") {
        par_init <- c(par_init, obj$c)
        lower_b <- c(lower_b, 0.0); upper_b <- c(upper_b, 0.45)
      }

      par_init <- c(par_init, obj$b)
      lower_b <- c(lower_b, rep(-8.0, length(obj$b)))
      upper_b <- c(upper_b, rep(8.0, length(obj$b)))

      # 3. Call Custom Solver
      res <- custom_solver(par_init, fn_neg_ll, lower_b, upper_b)

      # 4. Unpack Results
      p_new <- res$par
      ses <- get_se(res$hessian)
      idx <- 1

      if (obj$model %in% c("2PL","3PL","GPCM","GRM")) {
        obj$a <- p_new[idx]; obj$se_a <- ses[idx]; idx <- idx+1
      }
      if (obj$model == "3PL") {
        obj$c <- p_new[idx]; obj$se_c <- ses[idx]; idx <- idx+1
      }
      n_b <- length(obj$b)
      obj$b <- p_new[idx:(idx+n_b-1)]
      obj$se_b <- ses[idx:(idx+n_b-1)]

      item_list[[itm]] <- obj
    }

    diff <- abs(curr_ll - last_ll)
    if(con$verbose) message(sprintf("Cycle %d: LogLik = %.2f, Change = %.4f", cycle, curr_ll, diff))
    if (diff < con$conv_crit) converged <- TRUE
    last_ll <- curr_ll
  }

  # --- 5. Person Estimation (MLE) ---
  if(con$verbose) message("Estimating Person Parameters...")

  theta_est <- numeric(nrow(resp_mat))
  theta_se <- numeric(nrow(resp_mat))

  for (i in 1:nrow(resp_mat)) {
    obs_i <- resp_mat[i,]

    fn_person_neg_ll <- function(p_vec) {
      th <- p_vec[1]
      ll <- 0
      # Use item names directly to avoid index mismatches
      for (itm_name in items) {
        val <- obs_i[itm_name]
        if (!is.na(val)) {
          itm_obj <- item_list[[itm_name]]
          pr <- get_prob_matrix(th, itm_obj$a, itm_obj$b, itm_obj$c, itm_obj$model)
          # Ensure we don't exceed the number of columns in the prob matrix
          idx_val <- min(val + 1, ncol(pr))
          ll <- ll + log(pr[1, idx_val])
        }
      }
      return(-ll)
    }

    res <- custom_solver(c(0), fn_person_neg_ll, c(-6), c(6))
    theta_est[i] <- res$par[1]
    theta_se[i] <- get_se(res$hessian)[1]
  }

  # --- 6. Formatting ---
  out_rows <- list()
  for (itm in items) {
    obj <- item_list[[itm]]
    r <- list(
      item = itm, model = obj$model, status = ifelse(obj$type=="fixed","Fixed","Estimated"),
      discrimination = obj$a, discrimination_se = obj$se_a,
      guessing = if(obj$model=="3PL") obj$c else NA, guessing_se = obj$se_c
    )
    if (obj$model %in% c("Rasch","2PL","3PL")) {
      r$difficulty <- obj$b[1]; r$difficulty_se <- obj$se_b[1]
    } else {
      r$difficulty <- NA; r$difficulty_se <- NA
      for (k in seq_along(obj$b)) {
        r[[paste0("step_", k)]] <- obj$b[k]
        r[[paste0("step_", k, "_se")]] <- obj$se_b[k]
      }
    }
    out_rows[[itm]] <- r
  }

  all_cols <- unique(unlist(lapply(out_rows, names)))
  df_items <- data.frame(item = items)
  for (col in all_cols) if (col != "item") df_items[[col]] <- sapply(out_rows, function(x) if(col %in% names(x)) x[[col]] else NA)

  step_cols <- sort(grep("step_", names(df_items), value=TRUE))
  base_cols <- c("item", "model", "status", "discrimination", "discrimination_se",
                 "difficulty", "difficulty_se", "guessing", "guessing_se")
  df_items <- df_items[, c(base_cols, step_cols)[c(base_cols, step_cols) %in% names(df_items)]]

  df_fit <- data.frame(
    LogLik = last_ll,
    AIC = 2 * total_est_params - 2 * last_ll,
    BIC = total_est_params * log(nrow(resp_mat)) - 2 * last_ll,
    N_Est_Params = total_est_params, N_Students = nrow(resp_mat), Converged = converged
  )
  df_persons=data.frame(theta=theta_est,
                        theta_se=theta_se)
  message("---------------------------------------------------------")
  message("Process Complete. You can view results now.")

  return(list(
    item_params = df_items,
    person_params = df_persons,
    model_fit = df_fit
  ))
}
