pkgname <- "mutType"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('mutType')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("mutType")
### * mutType

flush(stderr()); flush(stdout())

### Name: mutType
### Title: Determine mutation type
### Aliases: mutType

### ** Examples


library(BSgenome.Hsapiens.UCSC.hg19)
Hs <- Hsapiens
sample <- system.file("extdata", "sample.vcf", package = "mutType")

mutType(sample, Hs, 7)




cleanEx()
nameEx("mutTypeTable")
### * mutTypeTable

flush(stderr()); flush(stdout())

### Name: mutTypeTable
### Title: Summarize Mutation Types
### Aliases: mutTypeTable

### ** Examples


library(BSgenome.Hsapiens.UCSC.hg19)
Hs <- Hsapiens
sample <- system.file("extdata", "sample.vcf", package = "mutType")

muts <- mutType(sample, Hs, 7, graphics = FALSE)
mutTypeTable(muts)





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
