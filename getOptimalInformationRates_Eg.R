#' @title Optimal timing for interim analyses in group sequential trials
#' @author Zhangyi He, Suzie Cro, Laurent Billot

#' OptimInterim (examples)

# call the getOptimalInformationRates function
source("./getOptimalInformationRates.R")

################################################################################

# Calculate optimal information rates in a four-stage group sequential design 
# (with an equal allocation and the O'Brien & Fleming type alpha- and beta-spending functions) 
# testing H0: mu1 - mu2 = 0 for an alternative H1: mu1 - mu2 = 0.2 
# with assumed standard deviation = 1; one-sided significance alpha = 0.025, power 1 - beta = 90%:
groups <- 2
alternative <- 0.2
stDev <- 1
allocationRatioPlanned <- 1
kMax <- 4
alpha <- 0.025
beta <- 0.1
sided <- 1
typeOfDesign <- "asOF"
typeBetaSpending <- "bsOF"

# get the optimal information rates in the group sequential design
informationRates <- getOptimalInformationRates(
  groups = groups,
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending)

# get the group sequential design with the optimal information rates
design <- getDesignGroupSequential(
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  informationRates = informationRates,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending)
summary(design)
print(design)
plot(design)
getDesignCharacteristics(design)

# get the sample sizes in the group sequential design with the optimal information rates
designPlan <- getSampleSizeMeans(
  design = design,
  groups = groups,
  alternative = alternative,
  stDev = stDev,
  allocationRatioPlanned = allocationRatioPlanned)
summary(designPlan)
print(designPlan)
plot(designPlan)

############################################################

# Calculate optimal information rates in a four-stage group sequential design 
# (with an equal allocation and O'Brien & Fleming type alpha- and beta-spending functions) 
# testing H0: pi1 - pi2 = 0 for an alternative H1: pi1 - pi2 = 0.2 
# with assumed pi2 = 0.2; one-sided significance alpha = 0.025, power 1 - beta = 90%:
groups <- 2
pi1 <- 0.4
pi2 <- 0.2
allocationRatioPlanned <- 1
kMax <- 4
alpha <- 0.025
beta <- 0.1
sided <- 1
typeOfDesign <- "asOF"
typeBetaSpending <- "bsOF"

# get the optimal information rates in the group sequential design
informationRates <- getOptimalInformationRatesRates(
  groups = groups,
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending)

# get the group sequential design with the optimal information rates
design <- getDesignGroupSequential(
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  informationRates = informationRates,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending)
summary(design)
print(design)
plot(design)
getDesignCharacteristics(design)

# get the sample sizes in the group sequential design with the optimal information rates
designPlan <- getSampleSizeRates(
  design = design,
  groups = groups,
  pi1 = pi1,
  pi2 = pi2,
  allocationRatioPlanned = allocationRatioPlanned)
summary(designPlan)
print(designPlan)
plot(designPlan)

################################################################################
