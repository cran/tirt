#' Simulate Testlet Response Theory Data (Vector Supported Version)
#'
#' @description
#' Simulate testlet responses data. Support both dichotomous and polytomous responses.
#' Provide an easy implementation with a few default settings.
#'
#' @param n_people Integer. Number of examinees.
#' @param item_structure List of lists defining item blocks.
#' @param theta Numeric vector (Optional). If provided, these exact ability values are used.
#' @param theta_mean Numeric. Mean of latent trait (used if theta is NULL).
#' @param theta_sd Numeric. SD of latent trait (used if theta is NULL).
#'
#' @return A list containing:
#'   \item{resp}{data.frame of responses (rows=people, cols=items)}
#'   \item{true_item_params}{data.frame of true item parameters}
#'   \item{true_person_params}{vector of true latent traits}
#' @examples
#'   # =========================================================================
#'   # Example 1: Complex Testlet Design
#'   # =========================================================================
#'   # Define the Testlet Blueprint
#'   trt_design <- list(
#'     # Testlet 1: Rasch Testlet Model (High dependence: var=0.8)
#'     list(model = "RaschT", n_items = 5, testlet_id = "Read_A",
#'          testlet_var = 0.8, b = c(-1, 1)),
#'
#'     # Testlet 2: 2PL Testlet Model (Default dependence: var=0.5)
#'     list(model = "2PLT", n_items = 5, testlet_id = "Read_B",
#'          a = c(0.7, 1.3)),
#'
#'     # Testlet 3: Graded Response Testlet (Polytomous, 4 categories)
#'     list(model = "GRT", n_items = 4, testlet_id = "Survey",
#'          categories = 4, testlet_var = 0.2)
#'   )
#'
#'   # Run Simulation
#'   trt_data <- sim_trt(n_people = 500, item_structure = trt_design)
#'
#'   # Inspect Results
#'   # 1. Responses
#'   head(trt_data$resp)
#'
#'   # 2. Item Parameters
#'   # (Notice 'testlet_loading' equals 'discrimination' for standard models)
#'   head(trt_data$true_item_params)
#'
#'   # 3. Person Parameters (Ability + Gamma for each testlet)
#'   head(trt_data$true_person_params)
#'
#'   # =========================================================================
#'   # Example 2: Manual Control (Theta, Gamma, and Parameters)
#'   # =========================================================================
#'
#'   # 1. Manual Theta (e.g., everyone has high ability)
#'   manual_theta <- rep(2.0, 100)
#'
#'   # 2. Manual Gamma (e.g., zero effect for T1)
#'   manual_gamma <- rep(0, 100)
#'
#'   # 3. Item Parameters: Exact Match vs Range Sampling
#'   custom_structure <- list(
#'     # Case A: Manual Gamma Vector
#'     list(model = "2PLT", n_items = 5, testlet_id = "T1",
#'          gamma_vector = manual_gamma),
#'
#'     # Case B: Exact Parameter Match (Length of 'a' equals n_items)
#'     list(model = "2PLT", n_items = 2, testlet_id = "T2",
#'          a = c(0.5, 2.5)),
#'
#'     # Case C: Range Sampling (Length of 'a' is 2, but n_items != 2)
#'     list(model = "2PLT", n_items = 5, testlet_id = "T3",
#'          a = c(0.5, 2.5))
#'   )
#'
#'   res_custom <- sim_trt(n_people = 100, theta = manual_theta,
#'                         item_structure = custom_structure)
#'
#'   # Verify Manual Theta
#'   print(mean(res_custom$true_person_params$ability)) # Should be 2.0
#'
#'   # Verify Manual Gamma (T1 should be 0)
#'   print(head(res_custom$true_person_params$testlet_T1))
#'
#'   # Verify Exact Match (T2 discrimination should be 0.5 and 2.5)
#'   print(res_custom$true_item_params[res_custom$true_item_params$testlet_id == "T2",
#'                                     "discrimination"])
#' @export
sim_trt <- function(n_people = 1000,
                    item_structure = list(),
                    theta = NULL,
                    theta_mean = 0,
                    theta_sd = 1) {

  # --- 1. Global Validation ---
  if (!is.numeric(n_people) || n_people < 1) stop("Error: 'n_people' must be a positive integer.")
  if (length(item_structure) == 0) stop("Error: 'item_structure' cannot be empty.")

  message("================================================================")
  message(sprintf("   STARTING TRT SIMULATION (N = %d)", n_people))
  message("================================================================")

  # --- 2. Initialize Person Parameters (Theta & Gamma) ---

  # A. Primary Trait (Theta) Logic
  if (!is.null(theta)) {
    # User provided vector
    if (length(theta) != n_people) {
      stop(sprintf("Error: Length of provided 'theta' vector (%d) does not match 'n_people' (%d).",
                   length(theta), n_people))
    }
    true_theta <- theta
    message(">> Ability (Theta): Using user-supplied vector.")
  } else {
    # Default Generation
    true_theta <- rnorm(n_people, mean = theta_mean, sd = theta_sd)
    message(sprintf(">> Ability (Theta): Generated from N(mean=%.2f, sd=%.2f).", theta_mean, theta_sd))
  }

  # B. Testlet Gammas
  testlet_ids <- c()
  testlet_configs <- list()

  # Pre-scan for configs
  for (i in seq_along(item_structure)) {
    blk <- item_structure[[i]]
    if (is.null(blk$testlet_id)) {
      tid <- paste0("T", i)
      item_structure[[i]]$testlet_id <- tid
      message(sprintf("Notice: Block %d missing 'testlet_id'. Assigned '%s'.", i, tid))
    } else {
      tid <- as.character(blk$testlet_id)
    }
    testlet_ids <- c(testlet_ids, tid)

    # Only configure if new
    if (is.null(testlet_configs[[tid]])) {
      if (!is.null(blk$gamma_vector)) {
        if (length(blk$gamma_vector) != n_people) stop(sprintf("Error (Block %d): 'gamma_vector' length mismatch.", i))
        testlet_configs[[tid]] <- list(type = "vector", val = blk$gamma_vector)
      } else if (!is.null(blk$testlet_var)) {
        testlet_configs[[tid]] <- list(type = "variance", val = blk$testlet_var)
      } else {
        testlet_configs[[tid]] <- list(type = "default", val = 0.5)
      }
    }
  }

  unique_tids <- unique(testlet_ids)
  gamma_matrix <- matrix(0, nrow = n_people, ncol = length(unique_tids))
  colnames(gamma_matrix) <- paste0("testlet_", unique_tids)

  message(">> Testlet Effects (Gamma):")
  for (j in seq_along(unique_tids)) {
    tid <- unique_tids[j]
    config <- testlet_configs[[tid]]

    if (config$type == "vector") {
      message(sprintf("   - Testlet '%s': Using user-supplied Gamma vector.", tid))
      gamma_matrix[, j] <- config$val
    } else if (config$type == "variance") {
      message(sprintf("   - Testlet '%s': Generated Gamma ~ N(0, %.2f) [User Var].", tid, config$val))
      gamma_matrix[, j] <- rnorm(n_people, mean = 0, sd = sqrt(config$val))
    } else {
      message(sprintf("   - Testlet '%s': Generated Gamma ~ N(0, 0.50) [Default].", tid))
      gamma_matrix[, j] <- rnorm(n_people, mean = 0, sd = sqrt(0.5))
    }
  }

  person_params <- data.frame(person_id = 1:n_people, ability = true_theta)
  person_params <- cbind(person_params, gamma_matrix)

  # --- 3. Item Simulation Loop ---
  response_list <- list()
  param_list <- list()
  current_item_idx <- 1

  valid_dich <- c("RaschT", "2PLT", "3PLT", "BiFT")
  valid_poly <- c("PCT", "GPCT", "GRT", "BiFT")
  all_valid <- c(valid_dich, valid_poly)

  for (i in seq_along(item_structure)) {
    block <- item_structure[[i]]
    if (is.null(block$model)) stop(sprintf("Block %d Error: Missing 'model'.", i))
    if (is.null(block$n_items)) stop(sprintf("Block %d Error: Missing 'n_items'.", i))

    model <- block$model
    n_items <- block$n_items
    tid <- as.character(block$testlet_id)

    if (!model %in% all_valid) stop(sprintf("Block %d Error: Invalid model '%s'.", i, model))

    # Categories Logic
    if (model %in% valid_dich && model != "BiFT") {
      cats <- 2
    } else {
      if (is.null(block$categories)) {
        if (model %in% valid_dich) cats <- 2 else {
          cats <- 3
          message(sprintf("   [Block %d] Categories not set. Defaulting to 3.", i))
        }
      } else {
        cats <- block$categories
        if (cats < 2) stop(sprintf("Block %d Error: Categories must be >= 2.", i))
      }
    }

    message(sprintf(">> Block %d: %d items (Model: %s, Testlet: %s)", i, n_items, model, tid))

    # --- SMART PARAMETER PARSER ---
    get_param <- function(name, default, n_req, friendly) {
      val <- block[[name]]

      # Case A: Not provided -> Default
      if (is.null(val)) {
        message(sprintf("      - %s: Default used (%.2f).", friendly, default))
        return(rep(default, n_req))
      }

      # Case B: Exact Vector Match (Highest Priority)
      # If user provides exactly n_items, we assume they want specific values per item
      if (length(val) == n_req) {
        if (n_req > 1) {
          message(sprintf("      - %s: Using user-supplied vector (Exact match).", friendly))
        } else {
          message(sprintf("      - %s: User fixed value (%.2f).", friendly, val[1]))
        }
        return(val)
      }

      # Case C: Range (Uniform Distribution)
      # Only if length is 2 AND n_req is NOT 2 (because Case B handles n=2)
      if (length(val) == 2 && val[1] < val[2]) {
        message(sprintf("      - %s: Sampling from U[%.2f, %.2f].", friendly, val[1], val[2]))
        return(runif(n_req, val[1], val[2]))
      }

      # Case D: Single Value (Fixed for all)
      if (length(val) == 1) {
        message(sprintf("      - %s: User fixed value (%.2f).", friendly, val))
        return(rep(val, n_req))
      }

      stop(sprintf("Block %d Error: Parameter '%s' has length %d but requires %d (or range definition).",
                   i, name, length(val), n_req))
    }

    # Generate Params
    if (model %in% c("RaschT", "PCT")) a <- rep(1, n_items) else a <- get_param("a", 1, n_items, "Discrimination a")
    if (model == "BiFT") s_load <- get_param("s", 1, n_items, "Testlet Loading s") else s_load <- a
    b <- get_param("b", 0, n_items, "Difficulty/Loc b")
    c_p <- if (model == "3PLT") get_param("c", 0, n_items, "Guessing c") else rep(0, n_items)

    block_params <- data.frame(
      item_id = paste0("item_", current_item_idx:(current_item_idx + n_items - 1)),
      model = model,
      testlet_id = tid,
      categories = cats,
      discrimination = a,
      testlet_loading = s_load,
      difficulty = b,
      guessing = c_p,
      stringsAsFactors = FALSE
    )

    gamma_col_name <- paste0("testlet_", tid)
    current_gamma <- person_params[[gamma_col_name]]
    block_resp <- matrix(NA, nrow = n_people, ncol = n_items)

    # Simulation Logic
    if (cats == 2) {
      for (j in 1:n_items) {
        z <- (a[j] * true_theta) + (s_load[j] * current_gamma) - (a[j] * b[j])
        prob <- c_p[j] + (1 - c_p[j]) / (1 + exp(-z))
        block_resp[, j] <- ifelse(runif(n_people) < prob, 1, 0)
      }
    } else {
      for (j in 1:n_items) {
        raw_steps <- sort(rnorm(cats - 1, mean = 0, sd = 0.5))
        thresholds <- b[j] + raw_steps
        for (k in 1:length(thresholds)) block_params[j, paste0("step_", k)] <- thresholds[k]

        if (model == "GRT" || (model == "BiFT" && cats > 2)) {
          prob_cum <- matrix(0, nrow = n_people, ncol = cats + 1)
          prob_cum[, 1] <- 1
          for (k in 1:(cats - 1)) {
            z <- (a[j] * true_theta) + (s_load[j] * current_gamma) - (a[j] * thresholds[k])
            prob_cum[, k + 1] <- 1 / (1 + exp(-z))
          }
          probs <- prob_cum[, 1:cats] - prob_cum[, 2:(cats+1)]
        } else {
          numerators <- matrix(0, nrow = n_people, ncol = cats)
          current_sum <- 0
          for (k in 1:(cats - 1)) {
            term <- (a[j] * true_theta) + (s_load[j] * current_gamma) - (a[j] * thresholds[k])
            current_sum <- current_sum + term
            numerators[, k + 1] <- current_sum
          }
          probs <- exp(numerators)
          probs <- probs / rowSums(probs)
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

  # --- 4. Final Assembly (Safe) ---
  message("================================================================")
  message("Constructing final dataframes...")

  all_resp <- as.data.frame(do.call(cbind, response_list))
  row.names(all_resp) <- NULL

  all_cols <- unique(unlist(lapply(param_list, names)))
  all_item_params <- do.call(rbind, lapply(param_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    if(length(missing) > 0) df[missing] <- NA
    return(df[, all_cols])
  }))
  row.names(all_item_params) <- NULL

  colnames(all_resp) <- all_item_params$item_id
  row.names(person_params) <- NULL

  message("Simulation Complete.")
  message(sprintf("Summary: %d items, %d examinees.", ncol(all_resp), n_people))
  message("================================================================")

  return(list(resp = all_resp, true_item_params = all_item_params, true_person_params = person_params))
}
