#' Differential Item Functioning (Mantel-Haenszel and Logistic Regression)
#'
#' @description
#' Detects Differential Item Functioning (DIF) for dichotomous items, that is,
#' items that behave differently for two groups of examinees (for example, a
#' reference and a focal group) after matching on ability. Two widely used
#' approaches are provided: the Mantel-Haenszel (1959) procedure, with the ETS
#' delta effect size and A/B/C flagging scheme (Holland & Thayer, 1988), and the
#' logistic regression approach (Swaminathan & Rogers, 1990), which separates
#' uniform and non-uniform DIF.
#'
#' @param data A data frame of dichotomous (0/1) item responses (rows = persons,
#'   columns = items).
#' @param group A vector of length equal to the number of rows in \code{data}
#'   giving the group membership of each person. It must have exactly two
#'   distinct non-missing values.
#' @param focal Optional. The value of \code{group} that identifies the focal
#'   group. If \code{NULL} (default), the second value (in sorted order) is used
#'   as the focal group and the first as the reference group.
#' @param match Optional numeric vector of matching scores (length equal to the
#'   number of rows in \code{data}). If \code{NULL} (default), the total number-
#'   correct score across all items is used.
#' @param purify Logical. If \code{TRUE}, the matching total score is computed
#'   as the rest score (total minus the studied item) separately for each item,
#'   which avoids contaminating the match with the item under study (default
#'   \code{FALSE}).
#'
#' @return A data frame with one row per item and the columns:
#' \itemize{
#'   \item \code{item}: the item name.
#'   \item \code{MH_chisq}, \code{MH_p}: the Mantel-Haenszel chi-square (with
#'     continuity correction) and its p-value.
#'   \item \code{MH_OR}: the Mantel-Haenszel common odds ratio.
#'   \item \code{MH_delta}: the ETS delta effect size, \eqn{-2.35 \ln(OR)}.
#'   \item \code{ETS_class}: the ETS DIF classification, \code{"A"} (negligible),
#'     \code{"B"} (moderate), or \code{"C"} (large).
#'   \item \code{LR_chisq}, \code{LR_p}: the 2-degree-of-freedom logistic
#'     regression test for combined (uniform + non-uniform) DIF.
#'   \item \code{LR_uniform_p}: p-value for uniform DIF (group main effect).
#'   \item \code{LR_nonuniform_p}: p-value for non-uniform DIF (group-by-match
#'     interaction).
#' }
#'
#' @references
#' Holland, P. W., & Thayer, D. T. (1988). Differential item performance and the
#' Mantel-Haenszel procedure. In H. Wainer & H. I. Braun (Eds.), \emph{Test
#' validity} (pp. 129-145). Erlbaum.
#'
#' Swaminathan, H., & Rogers, H. J. (1990). Detecting differential item
#' functioning using logistic regression procedures. \emph{Journal of
#' Educational Measurement, 27}(4), 361-370.
#'
#' @examples
#'   set.seed(123)
#'   sim <- sim_irt(n_people = 600,
#'                  item_structure = list(list(model = "2PL", n_items = 10)))
#'   resp <- sim$resp
#'
#'   # Create two groups and plant DIF in item 3 (harder for group "B")
#'   grp <- rep(c("A", "B"), each = 300)
#'   flip <- grp == "B" & resp[[3]] == 1
#'   resp[[3]][flip] <- rbinom(sum(flip), 1, 0.6)
#'
#'   dif_res <- dif(resp, group = grp)
#'   dif_res
#' @export
dif <- function(data,
                group,
                focal = NULL,
                match = NULL,
                purify = FALSE) {

  if (!is.data.frame(data) && !is.matrix(data)) stop("'data' must be a data frame or matrix.")
  resp <- as.matrix(data)
  N <- nrow(resp)
  J <- ncol(resp)
  item_names <- colnames(resp)
  if (is.null(item_names)) item_names <- paste0("Item_", seq_len(J))

  if (length(group) != N) stop("'group' must have the same length as the number of rows in 'data'.")

  # only dichotomous items are handled
  is_dich <- apply(resp, 2, function(x) all(x[!is.na(x)] %in% c(0, 1)))
  if (!all(is_dich)) {
    message(sprintf("Note: %d non-dichotomous item(s) detected and skipped. dif() supports dichotomous items.",
                    sum(!is_dich)))
  }

  grp_vals <- sort(unique(group[!is.na(group)]))
  if (length(grp_vals) != 2) stop("'group' must have exactly two distinct values.")
  if (is.null(focal)) focal <- grp_vals[2]
  if (!focal %in% grp_vals) stop("'focal' must be one of the values in 'group'.")
  ref <- setdiff(grp_vals, focal)

  is_focal <- group == focal
  is_ref   <- group == ref

  total_score <- rowSums(resp, na.rm = TRUE)
  if (is.null(match)) match <- total_score

  out <- vector("list", J)

  for (j in seq_len(J)) {
    y <- resp[, j]

    if (!is_dich[j]) {
      out[[j]] <- data.frame(item = item_names[j], MH_chisq = NA, MH_p = NA,
                             MH_OR = NA, MH_delta = NA, ETS_class = NA_character_,
                             LR_chisq = NA, LR_p = NA, LR_uniform_p = NA,
                             LR_nonuniform_p = NA, stringsAsFactors = FALSE)
      next
    }

    # matching variable (rest score if purify)
    mvar <- if (purify) total_score - y else match

    ok <- !is.na(y) & !is.na(mvar) & !is.na(group)

    # ---- Mantel-Haenszel over strata of the matching score ----
    strata <- split(which(ok), mvar[ok])
    num_or <- 0; den_or <- 0
    sum_a <- 0; sum_Ea <- 0; sum_Va <- 0

    for (idx in strata) {
      gy <- y[idx]; gg <- is_focal[idx]
      a1 <- sum(gy[!gg] == 1)   # reference correct
      b1 <- sum(gy[!gg] == 0)   # reference incorrect
      c1 <- sum(gy[gg]  == 1)   # focal correct
      d1 <- sum(gy[gg]  == 0)   # focal incorrect
      nk <- a1 + b1 + c1 + d1
      if (nk < 2) next
      num_or <- num_or + (a1 * d1) / nk
      den_or <- den_or + (b1 * c1) / nk

      nR <- a1 + b1; nF <- c1 + d1
      n1 <- a1 + c1; n0 <- b1 + d1
      sum_a  <- sum_a  + a1
      sum_Ea <- sum_Ea + nR * n1 / nk
      if (nk > 1) sum_Va <- sum_Va + (nR * nF * n1 * n0) / (nk^2 * (nk - 1))
    }

    mh_or <- if (den_or > 0) num_or / den_or else NA_real_
    mh_chisq <- if (sum_Va > 0) (abs(sum_a - sum_Ea) - 0.5)^2 / sum_Va else NA_real_
    mh_p <- if (is.na(mh_chisq)) NA_real_ else pchisq(mh_chisq, df = 1, lower.tail = FALSE)
    mh_delta <- if (!is.na(mh_or) && mh_or > 0) -2.35 * log(mh_or) else NA_real_

    # ETS A/B/C classification
    ets <- "A"
    if (!is.na(mh_delta) && !is.na(mh_p)) {
      sig <- mh_p < 0.05
      absd <- abs(mh_delta)
      if (sig && absd >= 1.5) ets <- "C"
      else if (sig && absd >= 1.0) ets <- "B"
      else ets <- "A"
    } else {
      ets <- NA_character_
    }

    # ---- Logistic regression DIF ----
    lr_chisq <- NA_real_; lr_p <- NA_real_; unif_p <- NA_real_; nonunif_p <- NA_real_
    dat_lr <- data.frame(y = y[ok], m = mvar[ok], g = factor(is_focal[ok]))
    if (length(unique(dat_lr$y)) == 2 && nrow(dat_lr) > 5) {
      m0 <- tryCatch(glm(y ~ m, data = dat_lr, family = binomial()), error = function(e) NULL)
      m1 <- tryCatch(glm(y ~ m + g, data = dat_lr, family = binomial()), error = function(e) NULL)
      m2 <- tryCatch(glm(y ~ m * g, data = dat_lr, family = binomial()), error = function(e) NULL)
      if (!is.null(m0) && !is.null(m2)) {
        lr_chisq <- as.numeric(m0$deviance - m2$deviance)
        lr_p <- pchisq(lr_chisq, df = 2, lower.tail = FALSE)
      }
      if (!is.null(m0) && !is.null(m1)) {
        unif_p <- pchisq(as.numeric(m0$deviance - m1$deviance), df = 1, lower.tail = FALSE)
      }
      if (!is.null(m1) && !is.null(m2)) {
        nonunif_p <- pchisq(as.numeric(m1$deviance - m2$deviance), df = 1, lower.tail = FALSE)
      }
    }

    out[[j]] <- data.frame(
      item = item_names[j],
      MH_chisq = round(mh_chisq, 3), MH_p = round(mh_p, 4),
      MH_OR = round(mh_or, 3), MH_delta = round(mh_delta, 3),
      ETS_class = ets,
      LR_chisq = round(lr_chisq, 3), LR_p = round(lr_p, 4),
      LR_uniform_p = round(unif_p, 4), LR_nonuniform_p = round(nonunif_p, 4),
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, out)
  row.names(res) <- NULL
  attr(res, "reference") <- ref
  attr(res, "focal") <- focal
  res
}
