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


#' generate a tibble of  MMRd TMB signature
#'
#' @param data4signeR data for signeR
#' @param data4sats data for signeR
#' @returns A tiblle
#' @import gtools
#' @import SATS
#' @import signeR
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import rtracklayer
#' @export
#' @examples
#' \dontrun{
#' mmrdSignature(data)
#' }
mmrdSignature <- function(data4signeR, data4sats) {
  opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                  target_regions_ocaplus, nsamples=nrow(data4signeR))
  H_hat <- EstimateSigActivity(V = data4sats, L = t(opp), W = W_TMB_MMR)
  SigExp <- CalculateSigExpectancy(L = t(opp), W = W_TMB_MMR, H = H_hat$H)
  SigExp <- round(SigExp, 2) |>
    t() |>
    as_tibble() |>
    mutate(sample =colnames(SigExp))
  return(SigExp)
}

#' generate a tibble of  Cosmic TMB signature
#'
#' @param data4signeR data for signeR
#' @param data4sats data for signeR
#' @returns A list including opp, H_hat, reconstruction profile, and SigExp
#' @import gtools
#' @import SATS
#' @import signeR
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import rtracklayer
#' @export
#' @examples
#' \dontrun{
#' cosmicSignature(data)
#' }

cosmicSignature <- function(data4signeR, data4sats) {
  opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                  target_regions_ocaplus, nsamples=nrow(data4signeR))
  H_hat <- EstimateSigActivity(V = data4sats, L = t(opp), W = W_TMB)
  SigExp <- CalculateSigExpectancy(L = t(opp), W = W_TMB, H = H_hat$H)
  SigExp <- round(SigExp, 2) |>
    t() |>
    as_tibble() |>
    mutate(sample =colnames(SigExp))
  reconProfile <- t(opp)*(as.matrix(W_TMB) %*% H_hat$H)
  return(list(opp, H_hat, reconProfile, SigExp))
}

