
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
#' @param context_length an odd integer representing the length of region
#' around the mutation
#' @param graphics default = TRUE; if TRUE, saves a barplot of mutation type
#' frequencies in a pdf file
#' @return data.frame object of mutations as "UP\[REF>ALT\]DOWN" and mutation
#' types as "[REF>ALT]" where "REF" is always either "C" or "T".
#' 
#' @examples
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' Hs <- Hsapiens
#' sample <- system.file("extdata", "sample.vcf", package = "mutType")
#' 
#' mutType(sample, Hs, 7)
#' 
#' @author Alessandro Giulivo\cr Politecnico di Milano\cr Maintainer:
#' Alessandro Giulivo\cr E-Mail: <alessandro.giulivo@@mail.polimi.it>
#' @references \url{https://github.com/AlessandroGiulivo/sciprog/tree/main/RProj}
#' 
#' @seealso \code{\link{mutTypeTable}}
#'
#' @export mutType
#' @importFrom VariantAnnotation readVcf
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom IRanges resize
#' @importFrom tidyr separate_rows
#' @importFrom dplyr filter
#' @import Biostrings
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @importFrom utils write.table

mutType <- function(VCFFile, refGenome, context_length, graphics = TRUE) {
  print('Reading file...', quote = FALSE)
  muts <- readVcf(VCFFile)
  muts <- rowRanges(muts)
  print('Computing regions...', quote = FALSE)
  muts <- resize(muts, width = context_length, fix="center")
  muts <- as.data.frame(muts, row.names = NULL)
  print('Checking chromosome names...', quote = FALSE)
  muts$seqnames <- apply(muts, 1,
                function (x) {ifelse(grepl('chr', x$seqnames),
                                    x$seqnames,
                                    as.character(paste0('chr', x$seqnames)))
                  }
                )
  muts$ALT <- apply(muts, 1, function(x) {ifelse(length(x$ALT) == 1,
                                                 as.character(unlist(x$ALT)),
                                                 paste(as.vector(unlist(x$ALT)),
                                                       collapse=','))})
  muts <- separate_rows(muts, ALT, sep=',', convert=TRUE)
  print('Filtering SNVs...', quote = FALSE)
  muts <- dplyr::filter(muts, nchar(REF) == 1 & nchar(ALT) == 1)
  muts <- as.data.frame(muts)
  print('Finding sequences...', quote = FALSE)
  muts$seq <- apply(muts, 1,
                    function (x)
                      refGenome[[x['seqnames']]][x['start']:x['end']])
  muts$up <- apply(muts, 1,
                   function (x)
                     substr(x$seq, 1, as.integer((context_length-1)/2)))
  muts$mut <- apply(muts, 1,
                    function (x)
                      paste0('[', x$REF, '>', x$ALT, ']'))
  muts$down <- apply(muts, 1,
                     function (x)
                       substr(x$seq,
                              as.integer((context_length+1)/2 + 1),
                              context_length))
  print('Checking redundancies...', quote = FALSE)
  muts$mutation <- apply(muts, 1,
                         function (x)
                           if(substr(x$mut, 2, 2)=='C' ||
                              substr(x$mut, 2, 2)=='T') {
                             paste0(x$up, x$mut, x$down)
                           } else {paste0(
                              reverseComplement(DNAString(x$down)),
                              '[',
                              complement(DNAString(x$REF)),
                              '>',
                              complement(DNAString(x$ALT)),
                              ']',
                              reverseComplement(DNAString(x$up))
                              )
                             }
                         )
  muts$mutType <- apply(muts, 1,
                        function (x)
                          {ifelse(substr(x$mut, 2, 2)=='C' |
                                    substr(x$mut, 2, 2)=='T',
                                  x$mut,
                                  paste0('[',
                                         complement(DNAString(x$REF)), '>',
                                         complement(DNAString(x$ALT)), ']'
                                         ))
                        }
  )
  muts$pos <- apply(muts, 1, function (x) (x$start + (context_length - 1) / 2))
  muts$seqnames <- apply(muts, 1, function(x) as.character(x$seqnames))
  print('Computing summary and plotting the results...', quote = FALSE)
  if(graphics==TRUE) {
    mutTypeTable(muts)
  }
  muts <- muts[, c('seqnames', 'pos', 'mutation', 'mutType')]
  print('Writing results file...', quote = FALSE)
  write.table(muts, 'mutationTypes.txt', sep='\t', col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  return(muts)

}
