#' Mixed-Format English Language Arts (ELA) Assessment Data (Form 1)
#'
#' A dataset containing binary and polytomous responses for demonstration.
#'
#' @format A data frame with 52417 rows and 47 columns:
#' \itemize{
#'   \item \code{ITEM1} - \code{ITEM30}: Binary responses (0 = Incorrect, 1 = Correct).
#'   \item \code{ITEM31} - \code{ITEM45}: Polytomous responses (scored 0-5).
#'   \item \code{THETA}: Latent ability estimates.
#'   \item \code{COVARIATE}: Person-level background variable.
#' }
#' @source Tang, C., Xiong, J., & Engelhard, G. (2025). Identification of writing
#' strategies in educational assessments with an unsupervised learning
#' measurement framework. \emph{Education Sciences}, 15(7), 912.
#' \doi{10.3390/educsci15070912}
#' @examples
#' data(ela1)
#' head(ela1)
"ela1"

#' Mixed-Format English Language Arts (ELA) Assessment Data (Form 2)
#'
#' A smaller dataset containing item responses.
#'
#' @format A data frame with columns representing item responses.
#' \itemize{
#'   \item \code{ITEM1} - \code{ITEM7}: Binary responses (0 = Incorrect, 1 = Correct).
#'   \item \code{ITEM8}: Polytomous response (scored 0-2).
#'   \item \code{ITEM9}: Polytomous response (scored 0-5).
#'   \item \code{ITEM10}: Polytomous response (scored 0-5).
#' }
#' @source Tang, C., Xiong, J., & Engelhard, G. (2025). Identification of writing
#' strategies in educational assessments with an unsupervised learning
#' measurement framework. \emph{Education Sciences}, 15(7), 912.
#' \doi{10.3390/educsci15070912}
#' @examples
#' data(ela2)
#' head(ela2)
"ela2"
