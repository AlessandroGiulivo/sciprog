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


mutType("data/inst/sample.vcf", HSapiens, 3)




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
