
#' Determine mutation type
#'
#' This function compares a list of mutations in VCF format to a reference genome
#' and returns the mutation type along with \emph context_length \emph upstream
#' and downstream bases
#'
#' @usage mutType(VCFFile, refGenome, context_length)
#' @param VCFFile a set of mutations in VCF Format
#' @param refGenome Full reference genome sequences as provided by UCSC in Biostrings objects
#' @param context_length an odd integer representing the length of region around the mutation
#' @return mutation type as "UP\[REF>ALT\]DOWN"
#' @author Alessandro Giulivo\cr Politecnico di Milano\cr Maintainer: Alessandro
#' Giulivo\cr E-Mail: <alessandro.giulivo@@mail.polimi.it>
#' 
#' @examples
#'
#' mutType("data/inst/sample.vcf", HSapiens, 3)
#'
#' @export
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
mutType <- function(VCFFile, refGenome, context_length) {
  muts <- readVcf(VCFFile)
  muts <- rowRanges(muts)
  return(muts)
}

