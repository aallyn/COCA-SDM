#### A collection of helper functions

## Check for and install any missing packages, from Vikram Baliga (http://www.vikram-baliga.com)
package_check <- function(packages) {
  lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}

