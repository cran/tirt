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

#' Large-Scale Mixed-Format English Language Arts (ELA) Assessment Data (Form 3)
#'
#' A long format dataset containing binary and polytomous responses for demonstration.
#'
#' @format A data frame with 2,434,185 observations and 5 variables:
#' \itemize{
#'   \item{STUDENTID}: Unique student identifier.
#'   \item{FORMID}: Unique test form identifier.
#'   \item{ITEMID}: Unique item identifier, where _T1 indicates trait 1 and _T2 indicates trait 2.
#'   \item{SEQ}: Position of the item in a form.
#'   \item{SCORE}: Dichotomously and polytomously scored item responses (0,1,2...).
#' }
#' @source Tang, C., Xiong, J., & Engelhard, G. (2025). Identification of writing
#' strategies in educational assessments with an unsupervised learning
#' measurement framework. \emph{Education Sciences}, 15(7), 912.
#' \doi{10.3390/educsci15070912}
#' @examples
#' data(ela3)
#' head(ela3)
"ela3"

#' Large-Scale Mixed-Format English Language Arts (ELA) Assessment Data Testmap (Form 3 Testmap)
#'
#' Test map for the long format ela3 response data
#'
#' @format A data frame with 328 rows and 5 variables:
#' \itemize{
#'   \item{FORMID}: Unique test form identifier.
#'   \item{SEQ}: Position of the item in a form.
#'   \item{ITEMID}: Unique item identifier, where _T1 indicates trait 1 and _T2 indicates trait 2.
#'   \item{TYPE}: Type of an item, where OP indicates operational items and FT indicates field-test items.
#'   \item{MAX_SCORE}: Max score of that item.
#' }
#' @source Tang, C., Xiong, J., & Engelhard, G. (2025). Identification of writing
#' strategies in educational assessments with an unsupervised learning
#' measurement framework. \emph{Education Sciences}, 15(7), 912.
#' \doi{10.3390/educsci15070912}
#' @examples
#' data(ela3_testmap)
#' head(ela3_testmap)
"ela3_testmap"


