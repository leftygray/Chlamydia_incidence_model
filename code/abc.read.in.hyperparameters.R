### R script to read in our hyperparameters for prior specification

### Supposes existence of a subdirectory called 'data'

### Required package
library(xlsx) # (available from online R archive: CRAN)

## Read in hyperparameters
hyperparameters <- read.xlsx("data/priors.xlsx",sheetIndex=1,colIndex=c(2:4),rowIndex=c(2:18),header=T) # my file: variable definitions as indicated by row names in .xlsx file (& by usage below, of course!)

prior.mvn.inf <- c(as.numeric(hyperparameters[9,1]),as.numeric(hyperparameters[10,1]),as.numeric(hyperparameters[11,1]),as.numeric(hyperparameters[12,1]))
prior.mvn.screening <- c(as.numeric(hyperparameters[13,1]),as.numeric(hyperparameters[14,1]),as.numeric(hyperparameters[15,1]),as.numeric(hyperparameters[16,1]))


