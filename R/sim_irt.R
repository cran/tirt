#' Simulate Item Response Theory Data
#'
#' @description
#' Simulate item responses data. Support both dichotomous and polytomous responses.
#' Provide an easy implementation with a few default settings.
#'
#' @param n_people Integer. Number of students.
#' @param item_structure List of lists defining item blocks.
#' @param theta Numeric vector (Optional). If provided, these exact ability values are used.
#' @param theta_mean Numeric. Mean of latent trait (used if theta is NULL).
#' @param theta_sd Numeric. SD of latent trait (used if theta is NULL).
#'
#' @return A list containing:
#'   \item{resp}{data.frame of responses (rows=people, cols=items)}
#'   \item{true_params}{data.frame of true item parameters}
#'   \item{theta}{vector of true latent traits}
#' @examples
#'   # 1. Define the Test Blueprint
#'   # We want:
#'   # - 10 items using 2PL (medium difficulty)
#'   # - 5 items using 3PL (difficult, with guessing)
#'   # - 5 items using GPCM (4-point Likert scale)
#'   # - 5 items using GRM (5-point Likert scale)
#'
#'   my_test_structure <- list(
#'     # Block 1: 2PL
#'     list(model = "2PL", n_items = 10, a = c(0.8, 1.2), b = c(-1, 1)),
#'
#'     # Block 2: 3PL (Harder items, b from 1 to 2.5, fixing guessing at 0.2)
#'     list(model = "3PL", n_items = 5, a = c(1.0, 1.5), b = c(1.0, 2.5), c = 0.2),
#'
#'     # Block 3: GPCM (Polytomous, 4 categories 0-3)
#'     list(model = "GPCM", n_items = 5, categories = 4, a = c(0.7, 1.3), b = c(-1, 1)),
#'
#'     # Block 4: GRM (Polytomous, 5 categories 0-4)
#'     list(model = "GRM", n_items = 5, categories = 5, a = c(1.0, 2.0))
#'   )
#'
#'   # 2. Run the Simulation
#'   # Define N and a specific Theta vector
#'   N <- 2000
#'   theta_vec <- rnorm(N, 0, 2)
#'
#'   sim_data <- sim_irt(
#'     n_people = N,
#'     theta = theta_vec,
#'     item_structure = my_test_structure
#'   )
#'
#'   # 3. Inspect the Output
#'   # The Response Matrix
#'   head(sim_data$resp)
#'
#'   # The True Parameters (Useful for recovery studies)
#'   # Note how it aligns a, b, and threshold parameters (step_1, step_2...)
#'   head(sim_data$true_params)
#' @export
sim_irt <- function(n_people = 1000,
                    item_structure = list(),
                    theta = NULL,
                    theta_mean = 0,
                    theta_sd = 1) {

  # --- INTERNAL HELPER: SANITIZE OUTPUTS ---
  # Strips hidden attributes to prevent RStudio Viewer crashes
  sanitize_output <- function(df_input) {
    data_list <- as.list(df_input)
    df_clean <- data.frame(data_list, stringsAsFactors = FALSE)
    row.names(df_clean) <- NULL
    return(df_clean)
  }

  # --- 1. Global Validation ---
  if (!is.numeric(n_people) || n_people < 1) {
    stop("Error: 'n_people' must be a positive integer (e.g., 500).")
  }

  message("----------------------------------------------------------------")
  message(sprintf("Starting Simulation for N = %d examinees...", n_people))

  # --- 2. Generate Theta ---
  if (!is.null(theta)) {
    if (length(theta) != n_people) {
      stop(sprintf("Error: Length of provided 'theta' vector (%d) does not match 'n_people' (%d).",
                   length(theta), n_people))
    }
    true_theta <- theta
    message(">> Ability (Theta): Using user-supplied vector.")
  } else {
    true_theta <- rnorm(n_people, mean = theta_mean, sd = theta_sd)
    message(sprintf(">> Ability (Theta): Generated from N(mean=%.2f, sd=%.2f).", theta_mean, theta_sd))
  }

  # Storage
  response_list <- list()
  param_list <- list()
  current_item_idx <- 1

  # Allowed Models Registry
  valid_dich <- c("Rasch", "2PL", "3PL")
  valid_poly <- c("PCM", "GPCM", "GRM")
  all_valid_models <- c(valid_dich, valid_poly)

  # --- 3. Iterate Blocks ---
  for (i in seq_along(item_structure)) {
    block <- item_structure[[i]]

    # Validation A: Missing Critical Keys
    if (is.null(block$model)) stop(sprintf("Block %d Error: You must specify 'model'.", i))
    if (is.null(block$n_items)) stop(sprintf("Block %d Error: You must specify 'n_items'.", i))

    model <- block$model
    n_items <- block$n_items

    # Validation B: Invalid Model Name
    if (!model %in% all_valid_models) {
      stop(sprintf(
        "Block %d Error: Invalid model '%s'.\nAllowed models: %s",
        i, model, paste(all_valid_models, collapse = ", ")
      ))
    }

    # Validation C: Invalid Item Count
    if (!is.numeric(n_items) || n_items < 1) {
      stop(sprintf("Block %d Error: 'n_items' must be a positive integer.", i))
    }

    message(sprintf(">> Block %d: %d items using %s", i, n_items, model))

    # --- Categories Logic ---
    if (model %in% valid_dich) {
      cats <- 2
      if (!is.null(block$categories) && block$categories != 2) {
        warning(sprintf("Block %d Warning: Model '%s' is dichotomous. Ignoring categories=%s.", i, model, block$categories))
      }
    } else {
      # Polytomous Models
      if (is.null(block$categories)) {
        cats <- 3
        message("   - Categories: Not specified. Defaulting to 3.")
      } else {
        cats <- block$categories
        if (!is.numeric(cats) || cats < 2) stop(sprintf("Block %d Error: Categories must be integer >= 2.", i))
        message(sprintf("   - Categories: User specified %d.", cats))
      }
    }

    # Init Block Storage
    block_resp <- matrix(NA, nrow = n_people, ncol = n_items)

    # Initialize Parameter Storage
    block_params <- data.frame(
      item_id = paste0("item_", current_item_idx:(current_item_idx + n_items - 1)),
      block = i,
      model = model,
      categories = cats,
      discrimination = NA,
      guessing = NA,
      stringsAsFactors = FALSE
    )

    # --- Helper: Smart Parameter Parser ---
    get_param <- function(param_name, default_val, n_req, friendly_name) {
      val <- block[[param_name]]

      # Case 1: Default
      if (is.null(val)) {
        message(sprintf("   - %s: Default values used (Fixed at %s).", friendly_name, default_val))
        return(rep(default_val, n_req))
      }

      # Case 2: Exact Vector Match (Highest Priority)
      # If provided vector length equals required items, use it exactly.
      if (length(val) == n_req) {
        if (n_req > 1) message(sprintf("   - %s: Using user-supplied vector (Exact match).", friendly_name))
        else message(sprintf("   - %s: User fixed value (%.2f).", friendly_name, val[1]))
        return(val)
      }

      # Case 3: Range (Uniform Distribution)
      # Only if length is 2 AND n_req is NOT 2 (because Case 2 handles n=2)
      if (length(val) == 2 && val[1] < val[2]) {
        message(sprintf("   - %s: Sampling from U[%.2f, %.2f].", friendly_name, val[1], val[2]))
        return(runif(n_req, val[1], val[2]))
      }

      # Case 4: Fixed Single Value
      if (length(val) == 1) {
        message(sprintf("   - %s: User fixed value (%.2f).", friendly_name, val))
        return(rep(val, n_req))
      }

      stop(sprintf("Block %d Error: Parameter '%s' length mismatch (expected %d or range).", i, param_name, n_req))
    }

    # --- Model Logic ---

    # A. Dichotomous
    if (model %in% valid_dich) {
      a <- if (model == "Rasch") rep(1, n_items) else get_param("a", 1, n_items, "Discrimination (a)")
      b <- get_param("b", 0, n_items, "Difficulty (b)")
      c_p <- if (model == "3PL") get_param("c", 0, n_items, "Guessing (c)") else rep(0, n_items)

      block_params$discrimination <- a
      block_params$difficulty <- b
      block_params$guessing <- c_p

      for (j in 1:n_items) {
        z <- a[j] * (true_theta - b[j])
        prob <- c_p[j] + (1 - c_p[j]) / (1 + exp(-z))
        block_resp[, j] <- ifelse(runif(n_people) < prob, 1, 0)
      }

      # B. Polytomous (PCM/GPCM)
    } else if (model %in% c("PCM", "GPCM")) {
      a <- if (model == "PCM") rep(1, n_items) else get_param("a", 1, n_items, "Discrimination (a)")
      loc <- get_param("b", 0, n_items, "Location (b)")

      block_params$discrimination <- a
      # Note: We do not store single difficulty here, we store steps below

      for (j in 1:n_items) {
        raw_steps <- sort(rnorm(cats - 1, 0, 0.5))
        thresholds <- loc[j] + raw_steps

        for (k in 1:length(thresholds)) {
          block_params[j, paste0("step_", k)] <- thresholds[k]
        }

        numerators <- matrix(0, nrow = n_people, ncol = cats)
        current_sum <- 0
        for (k in 1:(cats - 1)) {
          current_sum <- current_sum + (a[j] * (true_theta - thresholds[k]))
          numerators[, k + 1] <- current_sum
        }
        probs <- exp(numerators)
        probs <- probs / rowSums(probs)

        cum_probs <- t(apply(probs, 1, cumsum))
        rand_vals <- runif(n_people)
        block_resp[, j] <- apply(cum_probs >= rand_vals, 1, function(x) match(TRUE, x)) - 1
      }

      # C. GRM
    } else if (model == "GRM") {
      a <- get_param("a", 1, n_items, "Discrimination (a)")
      block_params$discrimination <- a
      message("   - Thresholds: Generated random ordered thresholds.")

      for (j in 1:n_items) {
        thresholds <- sort(runif(cats - 1, -2, 2))
        for (k in 1:length(thresholds)) {
          block_params[j, paste0("step_", k)] <- thresholds[k]
        }

        prob_cum <- matrix(0, nrow = n_people, ncol = cats + 1)
        prob_cum[, 1] <- 1
        for (k in 1:(cats - 1)) {
          prob_cum[, k + 1] <- 1 / (1 + exp(-(a[j] * (true_theta - thresholds[k]))))
        }
        probs <- prob_cum[, 1:cats] - prob_cum[, 2:(cats+1)]

        cum_probs_final <- t(apply(probs, 1, cumsum))
        rand_vals <- runif(n_people)
        block_resp[, j] <- apply(cum_probs_final >= rand_vals, 1, function(x) match(TRUE, x)) - 1
      }
    }

    response_list[[i]] <- block_resp
    param_list[[i]] <- block_params
    current_item_idx <- current_item_idx + n_items
  }

  # --- 4. Final Assembly (Sanitized) ---
  message("----------------------------------------------------------------")
  message("Constructing final data frames...")

  # A. Responses (Ensure safe conversion from matrix list)
  all_responses_mat <- do.call(cbind, response_list)
  all_responses <- as.data.frame(all_responses_mat, stringsAsFactors = FALSE)

  # B. Parameters (Handle missing columns)
  all_cols <- unique(unlist(lapply(param_list, names)))
  all_params <- do.call(rbind, lapply(param_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if(length(missing_cols) > 0) df[missing_cols] <- NA
    return(df[, all_cols])
  }))

  # Align Names
  colnames(all_responses) <- all_params$item_id

  # C. Sanitize (Strip Attributes)
  final_resp <- sanitize_output(all_responses)
  final_params <- sanitize_output(all_params)

  message("Simulation Complete.")
  message(sprintf("Summary: %d items, %d examinees.", ncol(final_resp), n_people))
  message("----------------------------------------------------------------")

  return(list(resp = final_resp,
              true_params = final_params,
              theta = true_theta))
}
