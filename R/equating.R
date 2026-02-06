#' Item Response Theory Equating / Linking
#'
#' @description
#' Conducts item response theory scale linking using Mean-Mean, Mean-Sigma, and Stocking-Lord methods.
#' Supports mixed formats of both dichotomous and polytomous models.
#' Automatically detects anchor items and validates model consistency.
#'
#' @param base_params Data frame of reference item parameters (Form X).
#' @param new_params Data frame of new item parameters to be transformed (Form Y).
#' @param person_params (Optional) Data frame of person parameters from Form Y.
#' @param methods Character vector. Options: "Mean-Mean", "Mean-Sigma", "Stocking-Lord".
#'   If NULL, defaults to all three.
#'
#' @return A list containing three data frames:
#'   \item{transformed_item_params}{New items transformed to Base scale (with SEs).}
#'   \item{transformed_person_params}{New persons transformed to Base scale (if provided).}
#'   \item{linking_constants}{The A (slope) and B (intercept) constants for each method.}
#' @examples
#'   # ===========================================================================
#'   # Example: Equating Form Y (New) to Form X (Base)
#'   # ===========================================================================
#'   set.seed(123)
#'
#'   # 1. Generate "True" Base Parameters (Form X)
#'   # ---------------------------------------------------------------------------
#'   # 10 Common Items (Anchors) + 10 Unique Items
#'   # 2PL and GRM mixed
#'
#'   gen_item_params <- function(n, type="2PL") {
#'     if(type=="2PL") {
#'       data.frame(
#'         item = paste0("Item_", 1:n),
#'         model = "2PL",
#'         a = round(runif(n, 0.8, 1.5), 2),
#'         b = round(rnorm(n, 0, 1), 2),
#'         stringsAsFactors = FALSE
#'       )
#'     } else {
#'       # GRM with 3 thresholds
#'       d <- t(apply(matrix(rnorm(n*3, 0, 0.5), n, 3), 1, sort))
#'       df <- data.frame(
#'         item = paste0("Poly_", 1:n),
#'         model = "GRM",
#'         a = round(runif(n, 0.8, 1.5), 2),
#'         stringsAsFactors = FALSE
#'       )
#'       df <- cbind(df, setNames(as.data.frame(d), paste0("step_", 1:3)))
#'       df
#'     }
#'   }
#'
#'   # Anchors
#'   anchor_2pl <- gen_item_params(5, "2PL")
#'   anchor_grm <- gen_item_params(3, "GRM")
#'   # Unique Form X
#'   unique_x <- gen_item_params(5, "2PL")
#'   unique_x$item <- paste0("X_", unique_x$item)
#'
#'   base_params <- dplyr::bind_rows(anchor_2pl, anchor_grm, unique_x)
#'
#'   # 2. Generate "New" Form Y Parameters (with Scale Shift)
#'   # ---------------------------------------------------------------------------
#'   # Scale Transformation: Theta_base = 1.2 * Theta_new + 0.5
#'   # True Constants: A = 1.2, B = 0.5
#'   TRUE_A <- 1.2
#'   TRUE_B <- 0.5
#'
#'   # Transform Anchor Parameters to "New" scale (Inverse Logic)
#'   # a_new = a_base * A
#'   # b_new = (b_base - B) / A
#'
#'   anchor_2pl_new <- anchor_2pl
#'   anchor_2pl_new$a <- anchor_2pl$a * TRUE_A
#'   anchor_2pl_new$b <- (anchor_2pl$b - TRUE_B) / TRUE_A
#'
#'   anchor_grm_new <- anchor_grm
#'   anchor_grm_new$a <- anchor_grm$a * TRUE_A
#'   step_cols <- grep("step_", names(anchor_grm_new))
#'   anchor_grm_new[, step_cols] <- (anchor_grm[, step_cols] - TRUE_B) / TRUE_A
#'
#'   # Unique Form Y
#'   unique_y <- gen_item_params(5, "2PL")
#'   unique_y$item <- paste0("Y_", unique_y$item)
#'
#'   new_params <- dplyr::bind_rows(anchor_2pl_new, anchor_grm_new, unique_y)
#'
#'   # 3. Create Dummy Person Parameters for Form Y
#'   # ---------------------------------------------------------------------------
#'   person_params <- data.frame(
#'     id = paste0("P", 1:50),
#'     theta = rnorm(50, 0, 1),
#'     theta_se = runif(50, 0.2, 0.5)
#'   )
#'
#'   # 4. Perform Equating
#'   # ---------------------------------------------------------------------------
#'   # We expect to recover A approx 1.2 and B approx 0.5
#'   results <- equate_irt(
#'     base_params = base_params,
#'     new_params = new_params,
#'     person_params = person_params,
#'     methods = c("Mean-Mean", "Stocking-Lord")
#'   )
#'
#'   # 5. Inspect Results
#'   # ---------------------------------------------------------------------------
#'   # Linking Constants
#'   print(results$linking_constants)
#'
#'   # Transformed Items (Form Y items on Form X scale)
#'   head(results$transformed_item_params)
#'
#'   # Transformed Persons
#'   head(results$transformed_person_params)
#' @export
equate_irt <- function(base_params, new_params, person_params = NULL,
                       methods = NULL) {

  # --- 1. Dependencies: Math Helpers ---

  # 1.1 Expected Score Calculator (for Stocking-Lord)
  get_expected_score <- function(theta, a, b_vec, c_param, model) {
    # Returns Expected Score (TCC contribution) for a single item at a vector of theta
    D <- 1.7

    # Binary
    if (model %in% c("Rasch", "2PL", "3PL")) {
      lin <- D * a * (theta - b_vec[1])
      p <- 1 / (1 + exp(-lin))
      if (model == "3PL") p <- c_param + (1 - c_param) * p
      return(p) # Expected score for binary is P(1)
    }

    # GPCM / PCM (Sum of k * Pk)
    if (model %in% c("GPCM", "PCM")) {
      n_cat <- length(b_vec) + 1
      numerators <- matrix(0, nrow = length(theta), ncol = n_cat)
      curr <- 0
      numerators[,1] <- 0
      for (k in 2:n_cat) {
        curr <- curr + a * (theta - b_vec[k-1])
        numerators[,k] <- curr
      }
      exps <- exp(numerators)
      probs <- exps / rowSums(exps)

      # Expected Score = Sum(k * Pk) where k=0..(m-1)
      # Pk is col k+1. Score is k.
      esc <- numeric(length(theta))
      for (k in 0:(n_cat-1)) {
        esc <- esc + k * probs[, k+1]
      }
      return(esc)
    }

    # GRM (Sum of Cumulative Probs)
    if (model == "GRM") {
      # ES = Sum_{k=1}^{m} P*_k
      esc <- numeric(length(theta))
      for (k in 1:length(b_vec)) {
        lin <- D * a * (theta - b_vec[k])
        p_star <- 1 / (1 + exp(-lin))
        esc <- esc + p_star
      }
      return(esc)
    }
    return(numeric(length(theta)))
  }

  # 1.2 Parameter Extractor
  extract_item_row <- function(row) {
    # Reads a DF row and outputs list(a, b_vec, c)
    mod <- row$model

    a <- if(!is.null(row$a) && !is.na(row$a)) row$a else 1.0
    if (mod %in% c("Rasch", "PCM")) a <- 1.0

    c_p <- if(!is.null(row$c) && !is.na(row$c)) row$c else 0.0
    if (!mod %in% c("3PL")) c_p <- 0.0

    # Find b/steps
    b_cols <- grep("^(b|d\\d+|step_\\d+)$", names(row), value=TRUE)
    b_cols <- b_cols[order(nchar(b_cols), b_cols)]
    b_vec <- as.numeric(row[b_cols])
    b_vec <- b_vec[!is.na(b_vec)]
    if(length(b_vec)==0 && "b" %in% names(row)) b_vec <- row$b

    return(list(a=a, b=b_vec, c=c_p, model=mod))
  }

  # --- 2. Input Validation & Setup ---

  message("---------------------------------------------------------")
  message("Starting IRT Equating / Linking")
  message("---------------------------------------------------------")

  # Default methods
  valid_methods <- c("Mean-Mean", "Mean-Sigma", "Stocking-Lord")
  if (is.null(methods)) methods <- valid_methods
  methods <- intersect(methods, valid_methods)

  if (length(methods) == 0) stop("No valid methods selected.")

  message(paste("Methods selected:", paste(methods, collapse=", ")))

  # Anchor Detection
  anchors <- intersect(base_params$item, new_params$item)
  n_anchors <- length(anchors)

  if (n_anchors == 0) stop("No common items (anchors) found between Base and New parameters.")
  if (n_anchors < 3) warning("Fewer than 3 anchor items detected. Linking may be unstable.")

  # Validate Models match
  base_anchors <- base_params[base_params$item %in% anchors, ]
  new_anchors <- new_params[new_params$item %in% anchors, ]

  # Sort to ensure alignment
  base_anchors <- base_anchors[order(base_anchors$item), ]
  new_anchors <- new_anchors[order(new_anchors$item), ]

  mismatches <- base_anchors$model != new_anchors$model
  if (any(mismatches)) {
    stop(paste("Model mismatch for anchor items:",
               paste(base_anchors$item[mismatches], collapse=", ")))
  }

  models_detected <- unique(c(base_params$model, new_params$model))
  message(paste("Models detected:", paste(models_detected, collapse=", ")))
  message(paste("Items detected: Base =", nrow(base_params), "| New =", nrow(new_params)))
  message(paste("Anchor items detected:", n_anchors))
  message(paste("Anchor IDs:", paste(anchors, collapse = ", ")))
  # ----------------------------
  if (!is.null(person_params)) message(paste("Persons to transform:", nrow(person_params)))

  # --- 3. Compute Linking Constants (A & B) ---

  results_list <- list() # Store A/B for each method

  # Pre-extract anchor stats for MM/MS
  # We need vectors of 'a' and 'b' (centroid)

  # Base stats
  ba_a <- numeric(n_anchors)
  ba_b <- numeric(n_anchors) # Centroid

  # New stats
  na_a <- numeric(n_anchors)
  na_b <- numeric(n_anchors)

  for (i in 1:n_anchors) {
    # Base
    p_b <- extract_item_row(base_anchors[i, ])
    ba_a[i] <- p_b$a
    ba_b[i] <- mean(p_b$b) # Mean of steps/difficulty

    # New
    p_n <- extract_item_row(new_anchors[i, ])
    na_a[i] <- p_n$a
    na_b[i] <- mean(p_n$b)
  }

  # --- Method 1: Mean-Mean ---
  if ("Mean-Mean" %in% methods) {
    message("Estimating constants: Mean-Mean...")
    # A = mean(a_base) / mean(a_new)
    # B = mean(b_base) - A * mean(b_new)

    # Note: For Rasch/PCM, a=1, so A=1.
    mu_a_base <- mean(ba_a)
    mu_a_new  <- mean(na_a)
    mu_b_base <- mean(ba_b)
    mu_b_new  <- mean(na_b)

    A_mm <- mu_a_base / mu_a_new
    B_mm <- mu_b_base - A_mm * mu_b_new

    results_list[["Mean-Mean"]] <- c(A = A_mm, B = B_mm)
  }

  # --- Method 2: Mean-Sigma ---
  if ("Mean-Sigma" %in% methods) {
    message("Estimating constants: Mean-Sigma...")
    # A = sd(b_base) / sd(b_new)
    # B = mean(b_base) - A * mean(b_new)

    sd_b_base <- sd(ba_b)
    sd_b_new  <- sd(na_b)

    # Safety for single item (SD=NA)
    if (is.na(sd_b_base) || is.na(sd_b_new) || sd_b_new == 0) {
      warning("Mean-Sigma requires >1 anchor item with variance. Skipping.")
    } else {
      A_ms <- sd_b_base / sd_b_new
      B_ms <- mean(ba_b) - A_ms * mean(na_b)
      results_list[["Mean-Sigma"]] <- c(A = A_ms, B = B_ms)
    }
  }

  # --- Method 3: Stocking-Lord ---
  if ("Stocking-Lord" %in% methods) {
    message("Estimating constants: Stocking-Lord (TCC Optimization)...")

    # Quadrature points for TCC integration
    nodes <- seq(-4, 4, length.out = 31)

    # Calculate Base TCC (Fixed)
    tcc_base <- numeric(length(nodes))
    for (i in 1:n_anchors) {
      p <- extract_item_row(base_anchors[i,])
      tcc_base <- tcc_base + get_expected_score(nodes, p$a, p$b, p$c, p$model)
    }

    # Optimization Function
    # Find A, B to transform New Params -> New TCC -> Match Base TCC
    loss_function <- function(par) {
      A <- par[1]
      B <- par[2]

      tcc_new_trans <- numeric(length(nodes))

      for (i in 1:n_anchors) {
        p <- extract_item_row(new_anchors[i,])

        # Transform parameters
        # a* = a / A
        # b* = A * b + B
        a_star <- p$a / A
        b_star <- A * p$b + B
        c_star <- p$c # Invariant

        tcc_new_trans <- tcc_new_trans + get_expected_score(nodes, a_star, b_star, c_star, p$model)
      }

      # Sum of squared differences
      return(sum((tcc_base - tcc_new_trans)^2))
    }

    # Optimize using L-BFGS-B (Allows bounds for A > 0)
    # Start A=1, B=0
    opt <- tryCatch({
      stats::optim(c(1, 0), loss_function, method = "L-BFGS-B",
                   lower = c(0.01, -10), upper = c(10, 10))
    }, error = function(e) {
      warning("Stocking-Lord optimization failed. Returning NA.")
      return(NULL)
    })

    if (!is.null(opt)) {
      results_list[["Stocking-Lord"]] <- c(A = opt$par[1], B = opt$par[2])
    }
  }

  # --- 4. Apply Transformations & Format Output ---

  message("Transforming parameters...")

  out_items <- data.frame()
  out_persons <- data.frame()
  out_consts <- data.frame()

  # Loop through calculated constants
  for (m in names(results_list)) {
    consts <- results_list[[m]]
    A <- consts["A"]
    B <- consts["B"]

    # 4.1 Transform Linking Constants DF
    out_consts <- rbind(out_consts, data.frame(
      Method = m, A = A, B = B, stringsAsFactors = FALSE
    ))

    # 4.2 Transform Items (New Form)
    for (i in 1:nrow(new_params)) {
      row <- new_params[i, ]
      p <- extract_item_row(row)

      # Transformation Logic
      # a* = a / A
      # b* = A * b + B
      # SE(a*) = SE(a) / A
      # SE(b*) = SE(b) * A

      # Transform Value
      new_a <- p$a / A
      new_b <- A * p$b + B
      new_c <- p$c

      # Transform SE (if exists)
      se_a <- if(!is.null(row$discrimination_se)) row$discrimination_se / A else NA
      se_c <- if(!is.null(row$guessing_se)) row$guessing_se  else NA #SE(c) remains invariant (unchanged).

      # Handle b SEs (need to match step columns)
      se_b_cols <- grep("(_se)$", names(row), value=TRUE)
      # Filter only difficulty/step SEs
      se_b_cols <- se_b_cols[!grepl("discrimination|guessing", se_b_cols)]
      se_b_vals <- as.numeric(row[se_b_cols])
      new_se_b <- se_b_vals * A

      # Construct Result Row
      res <- list(
        item = row$item,
        anchor_status = ifelse(row$item %in% anchors, "Anchor", "New"),
        method = m,
        model = row$model,
        discrimination = new_a,
        discrimination_se = se_a,
        guessing = new_c,
        guessing_se = se_c
      )

      # Add dynamic steps
      if (length(new_b) == 1 && row$model %in% c("Rasch", "2PL", "3PL")) {
        res$difficulty <- new_b[1]
        res$difficulty_se <- if(length(new_se_b)>0) new_se_b[1] else NA
      } else {
        res$difficulty <- NA
        for (k in 1:length(new_b)) {
          res[[paste0("step_", k)]] <- new_b[k]
          res[[paste0("step_", k, "_se")]] <- if(k <= length(new_se_b)) new_se_b[k] else NA
        }
      }

      # Bind to main DF (inefficient but safe for variable columns)
      # We convert list to single-row DF
      res_df <- as.data.frame(res, stringsAsFactors = FALSE)
      out_items <- dplyr::bind_rows(out_items, res_df) # Using dplyr bind_rows logic if avail, else base rbind
    }

    # 4.3 Transform Persons (if provided)
    if (!is.null(person_params)) {
      # theta* = A * theta + B
      # SE(theta*) = A * SE(theta)
      p_df <- person_params
      p_df$theta_trans <- A * p_df$theta + B
      p_df$theta_se_trans <- A * p_df$theta_se
      p_df$method <- m

      out_persons <- rbind(out_persons, p_df)
    }
  }

  # --- 5. Clean Up Output Columns ---

  # Reorder Item DF
  base_cols <- c("item", "anchor_status", "method", "model",
                 "discrimination", "discrimination_se",
                 "difficulty", "difficulty_se", "guessing", "guessing_se")
  step_cols <- sort(grep("step_", names(out_items), value=TRUE))
  final_cols <- c(base_cols, step_cols)
  final_cols <- intersect(final_cols, names(out_items))

  out_items <- out_items[, final_cols]
  row.names(out_items)<-NULL
  row.names(out_persons)<-NULL
  row.names(out_consts)<-NULL

  message("Equating completed.")
  message("---------------------------------------------------------")

  return(list(
    transformed_item_params = out_items,
    transformed_person_params = out_persons,
    linking_constants = out_consts
  ))
}
