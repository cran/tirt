#' Fixed Person Calibration with or without Covariate
#'
#' @description
#' Estimates item parameters (difficulty, discrimination) given fixed person parameters (theta),
#' with an optional person-level covariate. Supports Rasch and 2-Parameter Logistic models.
#'
#' @param df A data frame of item responses (0/1). Columns represent items, rows represent persons.
#' @param theta A numeric vector of person abilities (fixed parameters). Must match the number of rows in `df`.
#' @param model A character string specifying the model type. Options are "Rasch" or "2PL".
#' @param covariate An optional numeric vector representing a person-level covariate (e.g., time, group).
#'   Defaults to `NULL`.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item Item statistics (difficulty, standard errors, z-values, p-values).
#'   \item Discrimination parameters (for 2PL model).
#'   \item Global covariate effect (if `covariate` is provided).
#'   \item Classical item statistics (p-value, count, point-biserial correlation).
#'   \item Mean theta per item (average ability of persons answering the item).
#'   \item Infit and Outfit statistics (for Rasch model only).
#' }
#'
#' @importFrom dplyr select mutate group_by summarise ungroup full_join inner_join left_join rename rowwise across everything any_of pull
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_dbl
#' @importFrom stats glm binomial cor vcov pnorm setNames
#' @importFrom gtools mixedorder
#'
#' @examples
#'   # --- Example: With Selected Package Data ---
#'   data("ela1", package = "tirt")
#'
#'   # Subset data for a manageable example
#'   # Select the first 500 examinees and 30 item responses
#'   df_real <- ela1[1:500, 1:30]
#'
#'   # Extract pre-estimated latent traits and covariates
#'   fixed_theta <- ela1$THETA[1:500]
#'   fixed_cov <- ela1$COVARIATE[1:500]
#'
#'   # Estimate item parameters given fixed ability levels
#'   # fitting a 2-parameter logistic (2PL) model
#'   real_res <- fix_person(df = df_real,
#'                                    theta = fixed_theta,
#'                                    model = "2PL",
#'                                    covariate = fixed_cov)
#'   head(real_res)
#' \donttest{
#'   # --- Example: With Package Data ---
#'   data("ela1", package = "tirt")
#'
#'   # Select Item Responses (Cols 1-30)
#'   df_real <- ela1[, 1:30]
#'
#'   fixed_theta <- ela1$THETA
#'   fixed_cov <- ela1$COVARIATE
#'
#'   real_res <- fix_person(df = df_real,
#'                                    theta = fixed_theta,
#'                                    model = "2PL",
#'                                    covariate = fixed_cov)
#'   head(real_res)
#' }
#' @export
fix_person <- function(df, theta, model = c("Rasch", "2PL"), covariate = NULL) {

  # --- 1. Input Validation ---
  model <- match.arg(model)

  if (!is.data.frame(df)) stop("'df' must be a data frame.")
  if (length(theta) != nrow(df)) stop("Length of 'theta' must match the number of rows in 'df'.")
  if (!is.null(covariate) && length(covariate) != nrow(df)) stop("Length of 'covariate' must match the number of rows in 'df'.")

  # --- 2. Calculate Classical Statistics (Pre-cleaning) ---
  # Calculate item means (p-values) and counts
  df_stats <- df |>
    dplyr::summarise(dplyr::across(dplyr::everything(), list(
      pvalue = ~mean(.x, na.rm = TRUE),
      number = ~sum(!is.na(.x))
    ))) |>
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = c("item", ".value"),
      names_pattern = "(.+)_(.+)"
    )

  # Calculate point-biserial correlations with Theta
  # purrr::map_dbl works natively with |> as the first arg is .x (the data)
  item_cors <- df |>
    purrr::map_dbl(~ stats::cor(.x, theta, use = "complete.obs"))

  df_stats <- df_stats |>
    dplyr::mutate(correlation = item_cors[item])

  # Identify items to remove (deterministic 0 or 1)
  items_to_remove <- df_stats |>
    dplyr::filter(!is.na(pvalue), pvalue %in% c(0, 1)) |>
    dplyr::pull(item)

  if (length(items_to_remove) > 0) {
    warning(paste("The following items were removed due to p-value of 0 or 1:",
                  paste(items_to_remove, collapse = ", ")))
  }

  # --- 3. Data Preparation for GLM ---
  # Create a lightweight ID for joining
  df_aug <- df |>
    dplyr::mutate(
      .row_id = dplyr::row_number(),
      theta = theta
    )

  if (!is.null(covariate)) {
    df_aug$covariate <- covariate
  }

  # Pivot to long format and filter bad items
  df_long <- df_aug |>
    tidyr::pivot_longer(
      cols = -c(dplyr::any_of(c(".row_id", "theta", "covariate"))),
      names_to = "item",
      values_to = "response"
    ) |>
    dplyr::filter(!item %in% items_to_remove, !is.na(response))

  # Calculate Mean Theta (and Mean Covariate) per item based on available data
  item_means_agg <- df_long |>
    dplyr::group_by(item) |>
    dplyr::summarise(
      mean_theta = mean(theta, na.rm = TRUE),
      mean_covariate = if (!is.null(covariate)) mean(covariate, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  if (is.null(covariate)) item_means_agg$mean_covariate <- NULL


  # --- 4. Model Estimation ---

  # Construct formula
  if (model == "Rasch") {
    # Rasch: difficulty is intercept, slope is fixed at 1 via offset
    f_str <- "response ~ 0 + item + offset(theta)"
  } else {
    # 2PL: difficulty derived from intercept, slope is interaction
    f_str <- "response ~ 0 + item + item:theta"
  }

  if (!is.null(covariate)) {
    f_str <- paste(f_str, "+ covariate")
  }

  # Run GLM (base stats functions are not typically piped into)
  fit_glm <- stats::glm(
    as.formula(f_str),
    data = df_long,
    family = stats::binomial(link = "logit")
  )

  # --- 5. Parameter Extraction ---

  # Get raw coefficients
  coef_summ <- summary(fit_glm)$coefficients
  coef_df <- as.data.frame(coef_summ)
  coef_df$parameter <- rownames(coef_df)
  rownames(coef_df) <- NULL

  # Clean column names
  colnames(coef_df) <- c("Estimate", "Std_Error", "z_value", "Pr_z", "parameter")

  # Extract Covariate Effect (if exists)
  cov_row <- NULL
  if (!is.null(covariate)) {
    cov_row <- coef_df |>
      dplyr::filter(parameter == "covariate") |>
      dplyr::select(
        covariate = Estimate,
        covariate_se = Std_Error,
        covariate_zvalue = z_value,
        covariate_Pr_z = Pr_z
      )
  }

  # Process Item Parameters based on Model
  if (model == "Rasch") {

    # Filter item parameters
    result_items <- coef_df |>
      dplyr::filter(parameter != "covariate") |>
      dplyr::mutate(
        item = gsub("^item", "", parameter),
        # In Rasch GLM formulation: Beta = -Difficulty
        difficulty = -Estimate,
        difficulty_se = Std_Error,
        difficulty_zvalue = z_value,
        difficulty_Pr_z = Pr_z
      ) |>
      dplyr::select(item, difficulty, difficulty_se, difficulty_zvalue, difficulty_Pr_z)

  } else { # 2PL Model

    # Separate Intercepts (easiness) and Slopes (discrimination)
    # Pattern: Intercepts start with "item", Slopes contain ":theta"

    is_slope <- grepl(":theta", coef_df$parameter)
    is_intercept <- grepl("^item", coef_df$parameter) & !is_slope

    intercepts <- coef_df[is_intercept, ]
    slopes <- coef_df[is_slope, ]

    # Clean names for joining
    intercepts$item <- gsub("^item", "", intercepts$parameter)
    slopes$item <- gsub("^item|:theta", "", slopes$parameter)

    # Join
    item_table <- dplyr::inner_join(
      intercepts, slopes,
      by = "item",
      suffix = c("_int", "_slope")
    )

    # Calculate parameters and SEs
    vc <- stats::vcov(fit_glm)

    result_items <- item_table |>
      dplyr::rowwise() |>
      dplyr::mutate(
        discrimination = Estimate_slope,
        discrimination_se = Std_Error_slope,
        discrimination_zvalue = z_value_slope,
        discrimination_Pr_z = Pr_z_slope,

        difficulty = -Estimate_int / Estimate_slope,

        # Delta Method for SE(Difficulty)
        difficulty_se = {
          beta_name <- parameter_int
          a_name    <- parameter_slope

          beta   <- Estimate_int
          a      <- Estimate_slope

          # Gradient
          d_beta <- -1 / a
          d_a    <- beta / (a^2)
          grad   <- c(d_beta, d_a)

          # Covariance subset
          if (all(c(beta_name, a_name) %in% rownames(vc))) {
            cov_mat <- vc[c(beta_name, a_name), c(beta_name, a_name)]
            as.numeric(sqrt(t(grad) %*% cov_mat %*% grad))
          } else {
            NA_real_
          }
        },
        difficulty_zvalue = difficulty / difficulty_se,
        difficulty_Pr_z = 2 * (1 - stats::pnorm(abs(difficulty_zvalue)))
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        item,
        discrimination, discrimination_se, discrimination_zvalue, discrimination_Pr_z,
        difficulty, difficulty_se, difficulty_zvalue, difficulty_Pr_z
      )
  }

  # --- 6. Fit Statistics (Rasch Only) ---
  fit_stats_df <- NULL

  if (model == "Rasch") {

    # Prepare matrices for calculation
    O_mat <- as.matrix(df |> dplyr::select(dplyr::all_of(result_items$item)))

    diff_vec <- stats::setNames(result_items$difficulty, result_items$item)
    infit_vec <- numeric(length(diff_vec))
    outfit_vec <- numeric(length(diff_vec))
    item_names <- names(diff_vec)

    cov_effect <- 0
    if (!is.null(covariate) && !is.null(cov_row)) {
      cov_effect <- covariate * cov_row$covariate
    }

    # Loop over items
    for (j in seq_along(item_names)) {
      itm <- item_names[j]
      responses <- O_mat[, itm]

      lin_pred <- (theta - diff_vec[itm])
      if (!is.null(covariate)) {
        lin_pred <- lin_pred + cov_effect
      }

      E_j <- 1 / (1 + exp(-lin_pred))
      W_j <- E_j * (1 - E_j)

      valid_mask <- !is.na(responses)

      if (sum(valid_mask) > 0) {
        r_valid <- responses[valid_mask]
        e_valid <- E_j[valid_mask]
        w_valid <- W_j[valid_mask]

        resid_sq <- (r_valid - e_valid)^2

        infit_vec[j]  <- sum(w_valid * resid_sq) / sum(w_valid)
        outfit_vec[j] <- sum(resid_sq) / sum(w_valid)
      } else {
        infit_vec[j] <- NA
        outfit_vec[j] <- NA
      }
    }

    fit_stats_df <- data.frame(
      item = item_names,
      infit = infit_vec,
      outfit = outfit_vec,
      stringsAsFactors = FALSE
    )
  }

  # --- 7. Assembly and Output ---

  # Join Item Params + Stats + Theta Means
  final_df <- result_items |>
    dplyr::full_join(df_stats, by = "item") |>
    dplyr::full_join(item_means_agg, by = "item")

  # Join Covariate Global Params
  if (!is.null(cov_row)) {
    final_df <- cbind(final_df, cov_row)
  }

  # Join Fit Stats (if Rasch)
  if (!is.null(fit_stats_df)) {
    final_df <- final_df |>
      dplyr::left_join(fit_stats_df, by = "item")
  }

  # Sort by item naturally
  final_df <- final_df[gtools::mixedorder(final_df$item), ]

  # Define column order
  cols_base <- c("item")
  cols_2pl <- c("discrimination", "discrimination_se", "discrimination_zvalue", "discrimination_Pr_z")
  cols_diff <- c("difficulty", "difficulty_se", "difficulty_zvalue", "difficulty_Pr_z")
  cols_cov <- c("covariate", "covariate_se", "covariate_zvalue", "covariate_Pr_z")
  cols_stats <- c("pvalue", "number", "correlation", "mean_theta", "mean_covariate")
  cols_fit <- c("infit", "outfit")

  select_vec <- cols_base
  if (model == "2PL") select_vec <- c(select_vec, cols_2pl)
  select_vec <- c(select_vec, cols_diff)
  if (!is.null(covariate)) select_vec <- c(select_vec, cols_cov)
  select_vec <- c(select_vec, cols_stats)
  if (model == "Rasch") select_vec <- c(select_vec, cols_fit)

  # Final selection
  final_df <- final_df |>
    dplyr::select(dplyr::any_of(select_vec))

  rownames(final_df) <- NULL

  return(final_df)
}
