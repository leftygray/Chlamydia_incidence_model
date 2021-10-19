### R script to read in the key datasets for this project & to prepare them as R matrices

### Supposes existence of a subdirectory called 'data'

### Required package
# library(xlsx) # (available from online R archive: CRAN)

# Use a different library for reading in excel files so rJava is 
# not required
# load.library("readxl") 

## (1) ABS population estimates by (single-year) age & sex

# population.data <- read.xlsx("data/population.xlsx",sheetIndex=1,colIndex=c(1:202),rowIndex=c(1:13),header=F)
# population.data <- read_excel("data/population.xlsx")
# population.data <- read.csv("data/population-ng.csv")
population.data <- read.csv("data/population.csv")
population.data <- unname(as.matrix(population.data[,2:ncol(population.data)]))

nyears <- dim(population.data)[1] #+1 # Defines a global variable used throughout

m.population <- population.data[,1:101] #do.call(cbind,population.data[,1:101])
f.population <- population.data[,102:202] #do.call(cbind,population.data[,102:202])

# Population count in these single age groups used when computing incoming/outgoing numbers of infections across our 6 fixed age cohorts
m.14 <- m.population[,15] # column number = age - 1 (because of 0-1 age group)
m.19 <- m.population[,20]
m.24 <- m.population[,25]
m.29 <- m.population[,30]
m.34 <- m.population[,35]
f.14 <- f.population[,15]
f.19 <- f.population[,20]
f.24 <- f.population[,25]
f.29 <- f.population[,30]
f.34 <- f.population[,35]
m.14 <- c(m.14,m.14[12])
m.19 <- c(m.19,m.19[12])
m.24 <- c(m.24,m.24[12])
m.29 <- c(m.29,m.29[12])
m.34 <- c(m.34,m.34[12])
f.14 <- c(f.14,f.14[12])
f.19 <- c(f.19,f.19[12]) 
f.24 <- c(f.24,f.24[12]) 
f.29 <- c(f.29,f.29[12]) 
f.34 <- c(f.34,f.34[12])

# Population counts for our 6 fixed age cohorts
m.0.14 <- rowSums(m.population[,0:15])
m.15.19 <- rowSums(m.population[,16:20])
m.20.24 <- rowSums(m.population[,21:25])
m.25.29 <- rowSums(m.population[,26:30])
m.30.34 <- rowSums(m.population[,31:35])
m.35.plus <- rowSums(m.population[,36:101])
f.0.14 <- rowSums(f.population[,0:15])
f.15.19 <- rowSums(f.population[,16:20])
f.20.24 <- rowSums(f.population[,21:25])
f.25.29 <- rowSums(f.population[,26:30])
f.30.34 <- rowSums(f.population[,31:35])
f.35.plus <- rowSums(f.population[,36:101])
# m.0.14 <- c(m.0.14,m.0.14[12])
# m.15.19 <- c(m.15.19,m.15.19[12])
# m.20.24 <- c(m.20.24,m.20.24[12])
# m.25.29 <- c(m.25.29,m.25.29[12])
# m.30.34 <- c(m.30.34,m.30.34[12])
# m.35.plus <- c(m.35.plus,m.35.plus[12])
# f.0.14 <- c(f.0.14,f.0.14[12])
# f.15.19 <- c(f.15.19,f.15.19[12])
# f.20.24 <- c(f.20.24,f.20.24[12])
# f.25.29 <- c(f.25.29,f.25.29[12])
# f.30.34 <- c(f.30.34,f.30.34[12])
# f.35.plus <- c(f.35.plus,f.35.plus[12])

## (2) NNDSS notification and test counts for our fixed age & sex cohorts
# my file = years: 2001-2012 # 2006 & 2007: missing test counts
# nndss.data <- read.xlsx("data/notifications.xlsx",sheetIndex=1,colIndex=c(1:66),header=F) 
# nndss.data <- read_excel("data/notifications.xlsx", col_names = FALSE)
nndss.data <- read.csv("data/notifications.csv")
# nndss.data <- read.csv("data/notifications-ng.csv")

## Notification counts available for 0-14/15-19/20-24/25-29/30-34/35+
# notifications.m <- matrix(as.numeric(t(as.matrix(nndss.data[4:15,1:6]))),nrow=6) # this manipulation is needed because R reads this data initially as character strings (owing to the missing entries in the testing data)
notifications.m <- matrix(as.numeric(t(as.matrix(nndss.data[,2:7]))), nrow = 6)
# notifications.2013.m <- c(84,5081,11966,7890,3751,6068)
# notifications.new.m <- matrix(nrow=6,ncol=13)
# notifications.new.m[,1:12] <- notifications.m
# notifications.new.m[,13] <- notifications.2013.m
# notifications.m <- notifications.new.m
# notifications.f <- matrix(as.numeric(t(as.matrix(nndss.data[4:15,8:13]))),nrow=6)
notifications.f <- matrix(as.numeric(t(as.matrix(nndss.data[,9:14]))), nrow = 6)
# notifications.2013.f <- c(630,14561,18046,7865,3222,3320)
# notifications.new.f <- matrix(nrow=6,ncol=13)
# notifications.new.f[,1:12] <- notifications.f
# notifications.new.f[,13] <- notifications.2013.f
# notifications.f <- notifications.new.f

## Testing counts available for 0-14/15-24/25-34/35+
tests.data <- read.csv("data/tests.csv")
# tests.data <- read.csv("data/tests-ng.csv")

# tested.m <- matrix(0,nrow=4,ncol=12)
# for (i in 1:4) {for (j in 1:12) {tested.m[i,j] <- as.numeric(as.character(nndss.data[j+3,i+14]))}}
tested.m <- matrix(as.numeric(t(as.matrix(tests.data[,2:5]))), nrow = 4)
# tested.2013.m <- c(701,15515,14792,18478)+c(1102,54272,59283,63173)+c(770,24773,31153,34311)
# tested.new.m <- matrix(nrow=4,ncol=13)
# tested.new.m[,1:12] <- tested.m
# tested.new.m[,13] <- tested.2013.m
# tested.m <- tested.new.m
# tested.f <- matrix(0,nrow=4,ncol=12)
# for (i in 1:4) {for (j in 1:12) {tested.f[i,j] <- as.numeric(as.character(tests.data[j+3,i+19]))}}
tested.f <- matrix(as.numeric(t(as.matrix(tests.data[,7:10]))), nrow = 4)
# tested.2013.f <- c(1498,88914,71104,52895)+c(2848,161437,139050,103231)+c(1784,58845,51939,38837)
# tested.new.f <- matrix(nrow=4,ncol=13)
# tested.new.f[,1:12] <- tested.f
# tested.new.f[,13] <- tested.2013.f
# tested.f <- tested.new.f

