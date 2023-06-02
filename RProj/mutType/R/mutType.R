
#' Determine mutation type
#'
#' The mutType function compares a list of mutations in VCF format to a
#' reference genome and returns the mutation type along with upstream and
#' downstream bases for a total length = context_length
#'
#' @usage mutType(VCFFile, refGenome, context_length, graphics = TRUE)
#' @param VCFFile a set of mutations in VCF Format
#' @param refGenome Full reference genome sequences as provided by UCSC in
#' Biostrings objects
#' @param context_length an odd integer representing the length of region around
#' the mutation
#' @param graphics default = TRUE; if TRUE, saves a barplot of mutation type
#' frequencies in a pdf file
#' @return vector of mutation types as "UP\[REF>ALT\]DOWN"
#' 
#' @examples
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' Hs <- Hsapiens
#' sample <- system.file("extdata", "sample.vcf", package = "mutType")
#' 
#' mutType(sample, Hs, 7)
#' 
#' @author Alessandro Giulivo\cr Politecnico di Milano\cr Maintainer: Alessandro
#' Giulivo\cr E-Mail: <alessandro.giulivo@@mail.polimi.it>
#'
#' @export mutType
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges resize
#' @import Biostrings
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom utils write.table

mutType <- function(VCFFile, refGenome, context_length, graphics = TRUE) {
  muts <- readVcf(VCFFile)
  muts <- rowRanges(muts)
  muts <- resize(muts, width = context_length, fix="center")
  muts <- as.data.frame(muts)
  muts$seqnames <- lapply(muts$seqnames,
                          function (x) paste0('chr', x))
  muts$seq <- apply(muts, 1,
                    function (x) as.character(refGenome$chr19[x$start:x$end]))
  muts$up <- apply(muts, 1,
                   function (x) substr(x$seq, 0, as.integer((context_length-1)/2)))
  muts$mut <- apply(muts, 1,
                    function (x) paste0('[', x$REF, '>', x$ALT, ']'))
  muts$down <- apply(muts, 1,
                     function (x) substr(x$seq,
                                         as.integer((context_length+1)/2 + 1),
                                         context_length))
  muts$mutation <- apply(muts, 1,
                         function (x)
                           if(substr(x$mut, 2, 2)=='C' ||
                              substr(x$mut, 2, 2)=='T')
                             paste0(x$up, x$mut, x$down)
                            else paste0(
                              reverseComplement(DNAString(x$down)),
                              '[',
                              complement(DNAString(substr(x$mut, 2, 2))),
                              '>',
                              complement(DNAString(substr(x$mut, 4, 4))),
                              ']',
                              reverseComplement(DNAString(x$up))
                              )
                         )
  muts$mutType <- apply(muts, 1,
                        function (x) if(substr(x$mut, 2, 2)=='C' ||
                                        substr(x$mut, 2, 2)=='T') x$mut
                         else paste0('[',
                           complement(DNAString(substr(x$mut, 2, 2))), '>',
                           complement(DNAString(substr(x$mut, 4, 4))), ']'
                           )
  )
  muts$pos <- apply(muts, 1, function (x) (x$start + (context_length - 1) / 2))
  muts$ALT <- apply(muts, 1, function(x) substr(x$mut, 4, 4))
  muts$seqnames <- apply(muts, 1, function(x) as.character(x$seqnames))
  if(graphics==TRUE) {
    mutTypeTable(muts$mutType)
  }
  muts <- muts$mutType
  write.table(muts, 'mutationTypes.txt', sep='\n', col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  return(muts)

}
