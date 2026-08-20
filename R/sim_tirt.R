#' Simulate Mixed IRT and Testlet (TRT) Data
#'
#' @description
#' Simulates a data set that mixes ordinary independent item response theory
#' items with testlet (locally dependent) items in a single test form, the
#' simulation counterpart of \code{\link{irt_trt}}. Independent blocks depend
#' only on the primary trait (theta); testlet blocks additionally depend on a
#' testlet-specific nuisance effect (gamma). Both dichotomous and polytomous
#' formats are supported for the independent and the testlet parts, so the whole
#' family of models handled by \code{\link{irt_trt}} can be generated at once.
#'
#' @param n_people Integer. Number of examinees.
#' @param item_structure List of lists defining item blocks (see Details).
#' @param theta Numeric vector (Optional). If provided, these exact ability
#'   values are used.
#' @param theta_mean Numeric. Mean of the latent trait (used if \code{theta} is
#'   \code{NULL}).
#' @param theta_sd Numeric. SD of the latent trait (used if \code{theta} is
#'   \code{NULL}).
#'
#' @details
#' Each element of \code{item_structure} is a list describing one block. The
#' \code{model} name determines whether the block is independent or a testlet:
#' \itemize{
#'   \item Independent (theta only): \code{"Rasch"}, \code{"2PL"}, \code{"3PL"},
#'     \code{"PCM"}, \code{"GPCM"}, \code{"GRM"}.
#'   \item Testlet (theta + gamma): \code{"RaschT"}, \code{"2PLT"},
#'     \code{"3PLT"}, \code{"BiFT"}, \code{"PCT"}, \code{"GPCT"}, \code{"GRT"}.
#' }
#' Other recognized keys are \code{n_items}, \code{categories} (polytomous),
#' \code{a}, \code{b}, \code{c} (given as a fixed value or a \code{c(lo, hi)}
#' range to sample from), \code{s} (testlet loading for \code{"BiFT"}),
#' \code{testlet_id} (the testlet label, required for testlet blocks),
#' \code{testlet_var} (variance of the gamma effect), and \code{gamma_vector}
#' (a user-supplied gamma effect of length \code{n_people}).
#'
#' @return A list containing:
#'   \item{resp}{data.frame of responses (rows = people, cols = items).}
#'   \item{true_item_params}{data.frame of true item parameters, including a
#'     \code{model} column and a \code{testlet} column (\code{NA} for
#'     independent items). Its \code{item_id}, \code{model}, and \code{testlet}
#'     columns provide everything needed to build the \code{item_spec} argument
#'     of \code{\link{irt_trt}} (rename \code{item_id} to \code{item}).}
#'   \item{true_person_params}{data.frame of true person parameters: the primary
#'     ability plus one gamma column per testlet.}
#'
#' @seealso \code{\link{irt_trt}}, \code{\link{sim_trt}}, \code{\link{sim_irt}}
#'
#' @examples
#'   # A form with independent items and two testlets
#'   set.seed(2025)
#'   design <- list(
#'     list(model = "2PL",  n_items = 8),                       # independent
#'     list(model = "GRM",  n_items = 4, categories = 3),       # independent poly
#'     list(model = "2PLT", n_items = 4, testlet_id = "P1",
#'          testlet_var = 0.6),                                 # testlet
#'     list(model = "GPCT", n_items = 3, categories = 3,
#'          testlet_id = "P2", testlet_var = 0.5)               # testlet poly
#'   )
#'
#'   sim <- sim_tirt(n_people = 600, item_structure = design)
#'
#'   head(sim$resp)
#'   sim$true_item_params[, c("item_id", "model", "testlet")]
#'   head(sim$true_person_params)
#'
#'   \donttest{
#'   # Feed the true structure straight into the joint estimator
#'   spec <- data.frame(
#'     item = sim$true_item_params$item_id,
#'     model = sim$true_item_params$model,
#'     testlet = sim$true_item_params$testlet,
#'     stringsAsFactors = FALSE
#'   )
#'   fit <- irt_trt(sim$resp, spec, method = "EM",
#'                  control = list(max_iter = 15, verbose = FALSE))
#'   head(fit$item_params)
#'   }
#' @export
sim_tirt <- function(n_people = 1000,
                     item_structure = list(),
                     theta = NULL,
                     theta_mean = 0,
                     theta_sd = 1) {

  # --- 1. Validation ---
  if (!is.numeric(n_people) || n_people < 1) stop("Error: 'n_people' must be a positive integer.")
  if (length(item_structure) == 0) stop("Error: 'item_structure' cannot be empty.")

  message("================================================================")
  message(sprintf("   STARTING MIXED IRT/TRT SIMULATION (N = %d)", n_people))
  message("================================================================")

  # --- 2. Primary Trait (Theta) ---
  if (!is.null(theta)) {
    if (length(theta) != n_people) {
      stop(sprintf("Error: Length of provided 'theta' (%d) does not match 'n_people' (%d).",
                   length(theta), n_people))
    }
    true_theta <- theta
    message(">> Ability (Theta): Using user-supplied vector.")
  } else {
    true_theta <- rnorm(n_people, mean = theta_mean, sd = theta_sd)
    message(sprintf(">> Ability (Theta): Generated from N(mean=%.2f, sd=%.2f).", theta_mean, theta_sd))
  }

  independent_dich <- c("RASCH", "2PL", "3PL")
  independent_poly <- c("PCM", "GPCM", "GRM")
  testlet_dich <- c("RASCHT", "2PLT", "3PLT", "BIFT")
  testlet_poly <- c("PCT", "GPCT", "GPCMT", "GRT", "BIFT")

  is_testlet <- function(m) toupper(m) %in% c(testlet_dich, testlet_poly)

  # --- 3. Configure testlet gamma effects (only for testlet blocks) ---
  testlet_ids <- character(0)
  testlet_configs <- list()

  for (i in seq_along(item_structure)) {
    blk <- item_structure[[i]]
    if (is.null(blk$model)) stop(sprintf("Block %d Error: Missing 'model'.", i))
    if (!is_testlet(blk$model)) next

    if (is.null(blk$testlet_id)) {
      tid <- paste0("T", i)
      item_structure[[i]]$testlet_id <- tid
      message(sprintf("Notice: Testlet block %d missing 'testlet_id'. Assigned '%s'.", i, tid))
    } else {
      tid <- as.character(blk$testlet_id)
    }
    testlet_ids <- c(testlet_ids, tid)

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
  person_params <- data.frame(person_id = seq_len(n_people), ability = true_theta)

  if (length(unique_tids) > 0) {
    gamma_matrix <- matrix(0, n_people, length(unique_tids))
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
        gamma_matrix[, j] <- rnorm(n_people, 0, sqrt(config$val))
      } else {
        message(sprintf("   - Testlet '%s': Generated Gamma ~ N(0, 0.50) [Default].", tid))
        gamma_matrix[, j] <- rnorm(n_people, 0, sqrt(0.5))
      }
    }
    person_params <- cbind(person_params, gamma_matrix)
  } else {
    message(">> No testlet blocks detected: all items are independent.")
  }

  # --- 4. Item simulation loop ---
  response_list <- list()
  param_list <- list()
  current_item_idx <- 1

  for (i in seq_along(item_structure)) {
    block <- item_structure[[i]]
    if (is.null(block$n_items)) stop(sprintf("Block %d Error: Missing 'n_items'.", i))

    model <- block$model
    n_items <- block$n_items
    testlet_block <- is_testlet(model)
    tid <- if (testlet_block) as.character(block$testlet_id) else NA_character_
    mU <- toupper(model)

    # categories
    if (mU %in% c(independent_dich, testlet_dich) && mU != "BIFT") {
      cats <- 2
    } else {
      cats <- if (is.null(block$categories)) {
        if (mU %in% c(independent_dich, testlet_dich)) 2 else 3
      } else block$categories
      if (cats < 2) stop(sprintf("Block %d Error: 'categories' must be >= 2.", i))
    }

    message(sprintf(">> Block %d: %d items (Model: %s, %s)",
                    i, n_items, model,
                    if (testlet_block) paste0("Testlet: ", tid) else "Independent"))

    get_param <- function(name, default, n_req, friendly) {
      val <- block[[name]]
      if (is.null(val)) return(rep(default, n_req))
      if (length(val) == n_req) return(val)
      if (length(val) == 2 && val[1] < val[2]) return(runif(n_req, val[1], val[2]))
      if (length(val) == 1) return(rep(val, n_req))
      stop(sprintf("Block %d Error: Parameter '%s' length %d but requires %d (or a range).",
                   i, name, length(val), n_req))
    }

    a <- if (mU %in% c("RASCH", "RASCHT", "PCM", "PCT")) rep(1, n_items) else get_param("a", 1, n_items, "Discrimination a")
    s_load <- if (mU == "BIFT") get_param("s", 1, n_items, "Testlet Loading s") else a
    b <- get_param("b", 0, n_items, "Difficulty/Loc b")
    c_p <- if (mU %in% c("3PL", "3PLT")) get_param("c", 0, n_items, "Guessing c") else rep(0, n_items)

    gamma_vec <- if (testlet_block) person_params[[paste0("testlet_", tid)]] else rep(0, n_people)

    block_params <- data.frame(
      item_id = paste0("item_", current_item_idx:(current_item_idx + n_items - 1)),
      model = model,
      testlet = tid,
      categories = cats,
      discrimination = a,
      testlet_loading = if (testlet_block) s_load else rep(NA_real_, n_items),
      difficulty = b,
      guessing = c_p,
      stringsAsFactors = FALSE
    )

    block_resp <- matrix(NA, n_people, n_items)

    if (cats == 2) {
      for (j in seq_len(n_items)) {
        z <- (a[j] * true_theta) + (s_load[j] * gamma_vec) - (a[j] * b[j])
        prob <- c_p[j] + (1 - c_p[j]) / (1 + exp(-z))
        block_resp[, j] <- ifelse(runif(n_people) < prob, 1, 0)
      }
    } else {
      is_graded <- mU %in% c("GRM", "GRT") || (mU == "BIFT" && cats > 2)
      for (j in seq_len(n_items)) {
        raw_steps <- sort(rnorm(cats - 1, 0, 0.5))
        thresholds <- b[j] + raw_steps
        for (k in seq_along(thresholds)) block_params[j, paste0("step_", k)] <- thresholds[k]

        if (is_graded) {
          prob_cum <- matrix(0, n_people, cats + 1)
          prob_cum[, 1] <- 1
          for (k in seq_len(cats - 1)) {
            z <- (a[j] * true_theta) + (s_load[j] * gamma_vec) - (a[j] * thresholds[k])
            prob_cum[, k + 1] <- 1 / (1 + exp(-z))
          }
          probs <- prob_cum[, seq_len(cats)] - prob_cum[, seq_len(cats) + 1]
        } else {
          numer <- matrix(0, n_people, cats)
          running <- 0
          for (k in seq_len(cats - 1)) {
            running <- running + ((a[j] * true_theta) + (s_load[j] * gamma_vec) - (a[j] * thresholds[k]))
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

  # --- 5. Final assembly ---
  message("================================================================")
  message("Constructing final data frames...")

  all_resp <- as.data.frame(do.call(cbind, response_list))
  all_cols <- unique(unlist(lapply(param_list, names)))
  all_item_params <- do.call(rbind, lapply(param_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) df[missing_cols] <- NA
    df[, all_cols]
  }))
  row.names(all_resp) <- NULL
  row.names(all_item_params) <- NULL
  row.names(person_params) <- NULL
  colnames(all_resp) <- all_item_params$item_id

  message("Simulation Complete.")
  message(sprintf("Summary: %d items (%d testlets), %d examinees.",
                  ncol(all_resp), length(unique_tids), n_people))
  message("================================================================")

  list(resp = all_resp,
       true_item_params = all_item_params,
       true_person_params = person_params)
}
