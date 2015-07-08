# Package Creation Script
library(devtools)
# Sharing the package
rm(list = ls())
current.code <- as.package(".")
load_all(current.code)
document(current.code)
check(current.code) # I get an error but things still seem to work, so skipping for now
check_doc(current.code)

# if you want in your local environment to test
install(current.code)
library(CCPredict)
