## code to prepare `DATASET` dataset goes here

library(tidyverse)
library(SATS)
context <- read_delim("D:/rproject/mmrMsi/data/allVariants/M22-16155_LEEMS_c891_2023-07-25-10-53-29-832_2023-07-24_23-54-08-387_All/RESULTS/somatic_mutation_substitution_context_summary.txt")
data(SimData, package = "SATS")
SimData$V[1:4,1:4]
context4SATS <- context |>
  mutate(string = str_split(`#substitution_context`, ""))|>
  mutate(subs = map2_chr(string, substitution, function(x, y) str_c(x[1], "[", y, "]", x[3])))

rownames4SATS <- context4SATS$subs
#rownames(SimData$V) == colnames4SATS

mut <- read.table(system.file("extdata","21_breast_cancers.mutations.txt",
                              package="signeR"), header=TRUE, check.names=FALSE)

context
mut[1:3,1:3]

context4signeR <- context |>
  mutate(string = str_split(`#substitution_context`, ">"))|>
  mutate(subs = map2_chr(string, substitution, function(x, y) str_c(y, ":", x[1])))

colnames4signeR <- context4signeR$subs

#colnames(mut) == colnames4signeR

usethis::use_data(rownames4SATS, colnames4signeR, overwrite = TRUE)
usethis::use_data(target_regions)
