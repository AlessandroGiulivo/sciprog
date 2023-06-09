#' Summarize Mutation Types
#' 
#' The mutTypeTable function summarizes the mutations obtained from mutType
#' function into a counts table and produces a pdf image plotting the results.
#' It is automatically run when calling mutType function when
#' parameter graphics = T 
#'
#' @usage mutTypeTable(mutTypeResult)
#' @param mutTypeResult a data.frame as obtained from the mutType function
#' @return a table along with a histogram summarizing mutation types
#'
#' @examples
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' Hs <- Hsapiens
#' sample <- system.file("extdata", "sample.vcf", package = "mutType")
#' 
#' muts <- mutType(sample, Hs, 7, graphics = FALSE)
#' mutTypeTable(muts)
#'
#' @author Alessandro Giulivo\cr Politecnico di Milano\cr Maintainer:
#' Alessandro Giulivo\cr E-Mail: <alessandro.giulivo@@mail.polimi.it>
#' 
#' @references \url{https://github.com/AlessandroGiulivo/sciprog/tree/main/RProj}
#' 
#' @seealso \code{\link{mutType}}
#' 
#' @export mutTypeTable
#' @importFrom utils write.table
#' @import ggplot2 
mutTypeTable <- function(mutTypeResult) {
  tbl <- as.data.frame(table(mutTypeResult$mutType))
  colnames(tbl) <- c('Type', 'Frequency')
  p <- ggplot(tbl, aes(Type, Frequency, fill=Type)) +
    geom_col() +
    theme_minimal() +
    labs(title = 'Mutation Types Frequency', x='Mutation Type') +
    theme(legend.position = 'none') +
    scale_fill_viridis_d()
  ggsave("mutTypeFreqs.pdf", p)
  write.table(tbl, "mutationTypesTable.tsv", quote=FALSE,
              row.names=FALSE, sep='\t')
  print(p)
  return(tbl)
}
