#' @title Optimal timing for interim analyses in group sequential clinical trials
#' @author Zhangyi He, Laurent Billot, Suzie Cro

#' version 1.0

#' rpact plugin (examples)

setwd("~/Dropbox/Jeffery He/iResearch/Publications/2023/ZH2023-ClinTrial-StatAdv-GSD-OptimalInterim1")

# call the getOptimalInformationRates function
source("./Code/Code v1.0/getOptimalInformationRates_rpact_pluginR.R")

################################################################################

#' Get the optimal timing for interim analyses in the group sequential design for continuous endpoints

# O'Brien & Fleming spending

#
weights <- c(0, 0, 1, 0)
groups <- 2
thetaH0 <- 0
alternative <- 0.2
stDev <- 1
allocationRatioPlanned <- 1
kMax <- 4
alpha <- 0.05
beta <- 0.1
sided <- 1
typeOfDesign <- "asOF"
typeBetaSpending <- "bsOF"

# get the optimal timing for interim analyses in the group sequential design
system.time(res <- getOptimalInformationRates_ContEndpoint(
  weights = weights,
  groups = groups,
  thetaH0 = thetaH0,
  alternative = alternative,
  stDev = stDev,
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending))
informationRates <- res$informationRates
informationRates

# get the group sequential design with the optimal timing for interim analyses
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

# get the sample size in the group sequential design with the optimal timing for interim analyses
designPlan <- getSampleSizeMeans(
  design = design,
  groups = groups,
  alternative = alternative,
  stDev = stDev,
  allocationRatioPlanned = allocationRatioPlanned)
summary(designPlan)
print(designPlan)
plot(designPlan)

####################

#
weights <- c(0, 0, 1, 0)
groups <- 2
thetaH0 <- 0
alternative <- 0.2
stDev <- 1
allocationRatioPlanned <- NA
kMax <- 4
alpha <- 0.05
beta <- 0.1
sided <- 1
typeOfDesign <- "asOF"
typeBetaSpending <- "bsOF"

# get the optimal timing for interim analyses in the group sequential design
system.time(res <- getOptimalInformationRates_ContEndpoint(
  weights = weights,
  groups = groups,
  thetaH0 = thetaH0,
  alternative = alternative,
  stDev = stDev,
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending))
informationRates <- res$informationRates
informationRates
allocationRatioPlanned <- res$allocationRatioPlanned
allocationRatioPlanned

# get the group sequential design with the optimal timing for interim analyses
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

# get the sample size in the group sequential design with the optimal timing for interim analyses
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

#' Get the optimal timing for interim analyses in group sequential design for binary endpoints

# O'Brien & Fleming spending

#
weights <- c(0, 0, 1, 0)
groups <- 2
pi1 <- 0.4
pi2 <- 0.2
allocationRatioPlanned <- 1
kMax <- 4
alpha <- 0.05
beta <- 0.1
sided <- 1
typeOfDesign <- "asOF"
typeBetaSpending <- "bsOF"

# get the optimal timing for interim analyses in the group sequential design
system.time(res <- getOptimalInformationRates_BinEndpoint(
  weights = weights,
  groups = groups,
  pi1 = pi1,
  pi2 = pi2,
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending))
informationRates <- res$informationRates
informationRates

# get the group sequential design with the optimal timing for interim analyses
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

# get the sample size in the group sequential design with the optimal timing for interim analyses
designPlan <- getSampleSizeRates(
  design = design,
  groups = groups,
  pi1 = pi1,
  pi2 = pi2,
  allocationRatioPlanned = allocationRatioPlanned)
summary(designPlan)
print(designPlan)
plot(designPlan)

####################

#
weights <- c(0, 0, 1, 0)
groups <- 2
pi1 <- 0.4
pi2 <- 0.2
allocationRatioPlanned <- NA
kMax <- 4
alpha <- 0.05
beta <- 0.1
sided <- 1
typeOfDesign <- "asOF"
typeBetaSpending <- "bsOF"

# get the optimal timing for interim analyses in the group sequential design
system.time(res <- getOptimalInformationRates_BinEndpoint(
  weights = weights,
  groups = groups,
  pi1 = pi1,
  pi2 = pi2,
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending))
informationRates <- res$informationRates
informationRates
allocationRatioPlanned <- res$allocationRatioPlanned
allocationRatioPlanned

# get the group sequential design with the optimal timing for interim analyses
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

# get the sample size in the group sequential design with the optimal timing for interim analyses
designPlan <- getSampleSizeRates(
  design = design,
  groups = groups,
  pi1 = pi1,
  pi2 = pi2,
  allocationRatioPlanned = allocationRatioPlanned)
summary(designPlan)
print(designPlan)
plot(designPlan)

############################################################

#' Get the optimal timing for interim analyses in group sequential design for time-to-event endpoints



################################################################################
