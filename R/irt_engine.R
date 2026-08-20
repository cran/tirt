# =============================================================================
#  Internal computational engine (not exported)
# -----------------------------------------------------------------------------
#  These helpers translate the item-parameter data frames produced by the
#  estimation functions (binary_irt(), polytomous_irt(), mixed_irt(),
#  mirt_binary(), trt_binary(), trt_poly(), irt_trt(), equate_irt(),
#  fixed_item()) into a single normalized representation, and evaluate
#  category probabilities, item information, and expected item scores from it.
#
#  Keeping one shared engine guarantees that item_info(), test_info(),
#  person_fit(), score_table(), ld_stats(), item_fit(), tcc() and
#  reliability() all speak the same "language" as the estimators.
# =============================================================================

# Numerically safe logistic function.
.tirt_plogis <- function(x) {
  x <- pmin(pmax(x, -30), 30)
  1 / (1 + exp(-x))
}

# --- Identify the column holding item names -------------------------------
.tirt_item_column <- function(params) {
  nm <- names(params)
  for (cand in c("item", "item_id", "Item", "ITEM")) {
    if (cand %in% nm) return(as.character(params[[cand]]))
  }
  if (!is.null(rownames(params)) && !all(rownames(params) == seq_len(nrow(params)))) {
    return(rownames(params))
  }
  paste0("Item_", seq_len(nrow(params)))
}

# --- Locate a single-value column from a set of aliases -------------------
.tirt_pick <- function(params, aliases) {
  nm <- tolower(names(params))
  for (a in aliases) {
    hit <- which(nm == tolower(a))
    if (length(hit) > 0) return(params[[hit[1]]])
  }
  NULL
}

# --- Locate the ordered threshold/step columns ----------------------------
# Returns a character vector of column names (in category order), excluding
# any standard-error, z-value or p-value companion columns.
.tirt_threshold_cols <- function(params) {
  nm <- names(params)
  low <- tolower(nm)

  is_val <- grepl("^(thresh(old)?_[0-9]+|step_[0-9]+|step[/.]threshold_[0-9]+|d[0-9]+)$", low)
  is_se  <- grepl("(_se|_zvalue|_pr_z|_pr|_z)$", low)
  keep   <- is_val & !is_se

  cols <- nm[keep]
  if (length(cols) == 0) return(character(0))

  # order by the trailing integer
  ord <- as.integer(sub(".*?([0-9]+)$", "\\1", tolower(cols)))
  cols[order(ord)]
}

# --- Map any model label (incl. testlet variants) to a base family --------
.tirt_base_model <- function(model_label, n_cat) {
  if (is.null(model_label) || is.na(model_label) || model_label == "") {
    return(if (n_cat > 2) "GRM" else "2PL")
  }
  m <- toupper(as.character(model_label))
  m <- gsub("[^A-Z0-9]", "", m)
  # dichotomous families
  if (m %in% c("RASCH", "1PL", "RASCHT", "1PLT")) return("RASCH")
  if (m %in% c("2PL", "2PLT")) return("2PL")
  if (m %in% c("3PL", "3PLT")) return("3PL")
  # polytomous families
  if (m %in% c("GRM", "GRT")) return("GRM")
  if (m %in% c("GPCM", "GPCMT", "GPCT")) return("GPCM")
  if (m %in% c("PCM", "PCMT", "PCT")) return("PCM")
  # bifactor / unknown: decide from category count
  if (n_cat > 2) "GRM" else "2PL"
}

# =============================================================================
#  .tirt_parse_params()
#  Convert an item-parameter data frame into a list of per-item "specs".
#  Each spec is a list: list(model, a, b, c, n_cat) where
#    - b is the difficulty (binary) or an ordered vector of thresholds/steps.
#    - For GRM, b holds category thresholds; for GPCM/PCM, b holds the
#      cumulative category intercepts (the parameterization used internally
#      by polytomous_irt()/mixed_irt()).
# =============================================================================
.tirt_parse_params <- function(params, model = NULL, D = 1) {

  if (!is.data.frame(params)) {
    stop("'item_params' must be a data frame of item parameters ",
         "(for example, the $item_params element returned by binary_irt(), ",
         "polytomous_irt(), or mixed_irt()).")
  }

  item_names <- .tirt_item_column(params)
  J <- nrow(params)

  a_col <- .tirt_pick(params, c("discrimination", "a", "slope"))
  c_col <- .tirt_pick(params, c("guess", "guessing", "c"))
  b_col <- .tirt_pick(params, c("difficulty", "b", "location"))
  model_col <- .tirt_pick(params, c("model"))
  thr_cols  <- .tirt_threshold_cols(params)

  # User-supplied model override (single value or length-J vector)
  if (!is.null(model)) {
    if (length(model) == 1) model <- rep(model, J)
    if (length(model) != J) stop("'model' must have length 1 or match the number of items.")
  }

  specs <- vector("list", J)

  for (j in seq_len(J)) {
    a_j <- if (!is.null(a_col)) suppressWarnings(as.numeric(a_col[j])) else NA_real_
    c_j <- if (!is.null(c_col)) suppressWarnings(as.numeric(c_col[j])) else 0
    if (is.na(c_j)) c_j <- 0

    # gather non-missing threshold/step values for this item
    thr_j <- numeric(0)
    if (length(thr_cols) > 0) {
      raw <- suppressWarnings(as.numeric(unlist(params[j, thr_cols])))
      thr_j <- raw[!is.na(raw)]
    }

    if (length(thr_j) >= 1) {
      # --- polytomous item ---
      n_cat <- length(thr_j) + 1L
      lab <- if (!is.null(model))       model[j]
             else if (!is.null(model_col)) model_col[j]
             else NULL

      if (is.null(lab)) {
        # No explicit model label (e.g. polytomous_irt() output has no 'model'
        # column): infer the polytomous family from the column naming
        # convention -- 'step_' means GPCM/PCM, 'thresh_' means GRM.
        fam <- if (any(grepl("step", tolower(thr_cols)))) "GPCM" else "GRM"
      } else {
        fam <- .tirt_base_model(lab, n_cat)
        if (fam %in% c("RASCH", "2PL", "3PL")) {
          # label says binary but thresholds are present -> infer from names
          fam <- if (any(grepl("step", tolower(thr_cols)))) "GPCM" else "GRM"
        }
      }
      if (fam == "PCM") { fam <- "GPCM"; a_j <- 1 }
      if (is.na(a_j)) a_j <- 1

      specs[[j]] <- list(model = fam, a = a_j, b = thr_j, c = 0, n_cat = n_cat, D = D)

    } else {
      # --- dichotomous item ---
      b_j <- if (!is.null(b_col)) suppressWarnings(as.numeric(b_col[j])) else NA_real_
      lab <- if (!is.null(model))       model[j]
             else if (!is.null(model_col)) model_col[j]
             else NULL

      fam <- .tirt_base_model(lab, 2)
      if (fam %in% c("GRM", "GPCM", "PCM")) fam <- "2PL"  # label mismatch safety

      if (is.null(lab)) {
        # infer family from available columns
        if (!is.null(c_col) && !is.na(suppressWarnings(as.numeric(c_col[j])))) fam <- "3PL"
        else if (!is.null(a_col)) fam <- "2PL"
        else fam <- "RASCH"
      }
      if (fam == "RASCH") a_j <- 1
      if (is.na(a_j)) a_j <- 1
      if (fam != "3PL") c_j <- 0

      specs[[j]] <- list(model = fam, a = a_j, b = b_j, c = c_j, n_cat = 2L, D = D)
    }
  }

  list(item = item_names, spec = specs, n_items = J)
}

# =============================================================================
#  Category probabilities for one item spec at a vector of theta values.
#  Returns a matrix with length(theta) rows and n_cat columns.
# =============================================================================
.tirt_cat_probs <- function(spec, theta) {
  D <- spec$D
  a <- spec$a
  nq <- length(theta)

  if (spec$model %in% c("RASCH", "2PL", "3PL")) {
    b <- spec$b
    if (is.na(a) || is.na(b)) return(matrix(NA_real_, nq, 2))
    g <- .tirt_plogis(D * a * (theta - b))
    P <- spec$c + (1 - spec$c) * g
    return(cbind(1 - P, P))
  }

  if (spec$model == "GRM") {
    b <- spec$b
    K <- spec$n_cat
    Pstar <- matrix(0, nq, K + 1)
    Pstar[, 1] <- 1
    Pstar[, K + 1] <- 0
    bb <- sort(b)
    for (k in seq_len(K - 1)) {
      Pstar[, k + 1] <- .tirt_plogis(D * a * (theta - bb[k]))
    }
    probs <- Pstar[, seq_len(K), drop = FALSE] - Pstar[, seq_len(K) + 1, drop = FALSE]
    return(pmax(probs, 1e-12))
  }

  # GPCM / PCM (cumulative-intercept parameterization; step_0 = 0)
  K <- spec$n_cat
  steps <- c(0, spec$b)              # length K
  num <- matrix(0, nq, K)
  for (k in seq_len(K)) {
    num[, k] <- D * a * (k - 1) * theta - steps[k]
  }
  num <- num - apply(num, 1, max)   # stabilize
  ex <- exp(num)
  ex / rowSums(ex)
}

# =============================================================================
#  Fisher item information at a vector of theta values (returns a vector).
# =============================================================================
.tirt_item_info <- function(spec, theta) {
  D <- spec$D
  a <- spec$a

  if (is.na(a)) return(rep(NA_real_, length(theta)))

  if (spec$model %in% c("RASCH", "2PL")) {
    b <- spec$b
    if (is.na(b)) return(rep(NA_real_, length(theta)))
    g <- .tirt_plogis(D * a * (theta - b))
    return((D * a)^2 * g * (1 - g))
  }

  if (spec$model == "3PL") {
    b <- spec$b
    if (is.na(b)) return(rep(NA_real_, length(theta)))
    cc <- spec$c
    g <- .tirt_plogis(D * a * (theta - b))
    P <- cc + (1 - cc) * g
    P <- pmin(pmax(P, 1e-10), 1 - 1e-10)
    # Standard Birnbaum 3PL information; reduces to 2PL when c = 0.
    return((D * a)^2 * (1 - cc) * (1 - g) * g^2 / P)
  }

  if (spec$model == "GRM") {
    b <- sort(spec$b)
    K <- spec$n_cat
    nq <- length(theta)
    Pstar <- matrix(0, nq, K + 1)
    Pstar[, 1] <- 1
    Pstar[, K + 1] <- 0
    for (k in seq_len(K - 1)) Pstar[, k + 1] <- .tirt_plogis(D * a * (theta - b[k]))
    dPstar <- Pstar * (1 - Pstar) * (D * a)   # derivative of each cumulative prob
    dPstar[, 1] <- 0
    dPstar[, K + 1] <- 0
    probs <- pmax(Pstar[, seq_len(K), drop = FALSE] - Pstar[, seq_len(K) + 1, drop = FALSE], 1e-12)
    dprob <- dPstar[, seq_len(K), drop = FALSE] - dPstar[, seq_len(K) + 1, drop = FALSE]
    return(rowSums(dprob^2 / probs))
  }

  # GPCM / PCM : information = (D a)^2 * Var(score | theta)
  probs <- .tirt_cat_probs(spec, theta)
  scores <- 0:(spec$n_cat - 1)
  ex <- as.vector(probs %*% scores)
  ex2 <- as.vector(probs %*% (scores^2))
  (D * a)^2 * (ex2 - ex^2)
}

# =============================================================================
#  Expected item score at a vector of theta values (returns a vector).
# =============================================================================
.tirt_expected_score <- function(spec, theta) {
  probs <- .tirt_cat_probs(spec, theta)
  scores <- 0:(spec$n_cat - 1)
  as.vector(probs %*% scores)
}

# =============================================================================
#  Maximum possible score for an item (n_cat - 1).
# =============================================================================
.tirt_max_score <- function(spec) spec$n_cat - 1L

# =============================================================================
#  Normalize a supplied theta argument (vector, or a person-parameter data
#  frame with an 'ability'/'theta' column) into a plain numeric vector.
# =============================================================================
.tirt_extract_theta <- function(theta) {
  if (is.data.frame(theta)) {
    for (cand in c("ability", "theta", "Theta", "EAP")) {
      if (cand %in% names(theta)) return(as.numeric(theta[[cand]]))
    }
    stop("Could not find an 'ability' or 'theta' column in the supplied person parameters.")
  }
  as.numeric(theta)
}
