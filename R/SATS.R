#' generate a context matrix for SATS
#'
#' @param context context R object getting from OCA Plus file.
#' @returns A matirx.
#' @import dplyr
#' @import purrr
#' @import stringr
#' @export
#' @examples
#' \dontrun{
#' subsCount(context)
#' }
context4sats <- function(context) {
  subsCount <- context |>
    dplyr::select(counts)
  subsCount <- as.matrix(subsCount)
  rownames(subsCount) <- rownames4SATS

  return(subsCount)
}

#' generate a context matrix for SATS
#'
#' @param context context R object getting from OCA Plus file.
#' @returns A matirx.
#' @import dplyr
#' @import purrr
#' @import stringr
#' @export
#' @examples
#' \dontrun{
#' subsCount(context)
#' }
context4signeR <- function(context) {
  subsCount <- context |>
    dplyr::select(counts) |>
    t()
  subsCount <- as.matrix(subsCount, nrow = 1)
  colnames(subsCount) <- colnames4signeR
  return(subsCount)
}
