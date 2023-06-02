## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(mutType)

## ----mutType_example----------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
Hs <- Hsapiens

sample <- system.file("extdata", "sample.vcf", package = "mutType")

mutType(sample, Hs, 7, graphics = FALSE)

## ----mutTypeTable_example-----------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
Hs <- Hsapiens
sample <- system.file("extdata", "sample.vcf", package = "mutType")

muts <- mutType(sample, Hs, 7, graphics = FALSE)
mutTypeTable(muts)

