#' Fix a filename consistantly
#'
#' @param filename context R object getting from OCA Plus file.
#' @returns A file name characters
#' @import dplyr
#' @import purrr
#' @import stringr
#' @export
#' @examples
#' \dontrun{
#' renameFile(filename)
#' }
renameFile <- function(filename){
  id <- filename|>str_extract("[A-Z]{1,2}[0-9]{2}-[0-9]{1,6}")
  #id <- filenameSplits[1]
  idSplits <- str_split(id, "-") |> unlist()

  if(length(idSplits) != 2)
    stop("please check id format (M**-numbers)")
  if(str_length(idSplits[2]) > 6)
    stop("please check id format. The id number digit should be less than 6")

  lengthZero <- 6 - str_length(idSplits[2])

  zeros <- str_dup("0", lengthZero)
  newname <- str_c(idSplits[1], "-", zeros, idSplits[2], ".vcf")

  if(!grepl("M[0-9]{2}-[0-9]{6}\\.vcf", newname))
    stop("Check file name")

  return(newname)
}

#' generate a file list
#'
#' @param filePath context R object getting from OCA Plus file.
#' @returns a character vector
#' @export
#' @examples
#' \dontrun{
#' filesContext(here::here("path/to/contextfiles"))
#' }
filesContext <- function(filePath) {
  list.files(path = filePath, pattern = "somatic_mutation_substitution_context_summary", all.files = FALSE,
                           full.names = TRUE, recursive = TRUE,
                           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
}

#' generate a context matrix for SATS
#'
#' @param fileList file list generated from fileContext.
#' @returns a character vector
#' @import purrr
#' @export
#' @examples
#' \dontrun{
#' extractID(fileList)
#' }
extractID <- function(fileList){
  gsub(".vcf", "", fileList |> purrr::map(renameFile))
}
