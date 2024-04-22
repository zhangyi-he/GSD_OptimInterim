#' OptimInterim

#' Examples

# call the getOptimalInformationRates function
source("./OptimInterim.R")

################################################################################

# Calculate optimal information rates in a four-stage group sequential design
# (with an equal allocation and Haybittle-Peto boundaries)
# testing H0: mu1 - mu2 = 0 for an alternative H1: mu1 - mu2 = 0.2 with assumed standard deviation = 1;
# two-sided significance alpha = 0.05, power 1 - beta = 90%:
alternative <- 0.2
stDev <- 1
allocationRatioPlanned <- 1
kMax <- 4
alpha <- 0.05
beta <- 0.1
sided <- 2
typeOfDesign <- "HP"

# get the optimal information rates in the group sequential design
informationRates <- getOptimalInformationRates(
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign)
print(informationRates)

# get the group sequential design with the optimal information rates
design <- getDesignGroupSequential(
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  informationRates = informationRates,
  typeOfDesign = typeOfDesign)
summary(design)
print(design)
plot(design)
getDesignCharacteristics(design)

# get the sample sizes in the group sequential design with the optimal information rates
designPlan <- getSampleSizeMeans(
  design = design,
  alternative = alternative,
  stDev = stDev,
  allocationRatioPlanned = allocationRatioPlanned)
summary(designPlan)
print(designPlan)
plot(designPlan)

############################################################

# Calculate optimal information rates in a three-stage group sequential design
# (with an equal allocation and O'Brien & Fleming type alpha-spending boundaries)
# testing H0: pi1 - pi2 = 0 for an alternative H1: pi1 - pi2 = 0.2 with assumed pi2 = 0.2;
# one-sided significance alpha = 0.025, power 1 - beta = 80%:
pi1 <- 0.4
pi2 <- 0.2
allocationRatioPlanned <- 1
kMax <- 3
alpha <- 0.025
beta <- 0.2
sided <- 1
typeOfDesign <- "asOF"

# get the optimal information rates in the group sequential design
informationRates <- getOptimalInformationRates(
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign)
print(informationRates)

# get the group sequential design with the optimal information rates
design <- getDesignGroupSequential(
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  informationRates = informationRates,
  typeOfDesign = typeOfDesign)
summary(design)
print(design)
plot(design)
getDesignCharacteristics(design)

# get the sample sizes in the group sequential design with the optimal information rates
designPlan <- getSampleSizeRates(
  design = design,
  pi1 = pi1,
  pi2 = pi2,
  allocationRatioPlanned = allocationRatioPlanned)
summary(designPlan)
print(designPlan)
plot(designPlan)

############################################################

# Calculate optimal information rates in a four-stage group sequential design
# (with an equal allocation and O'Brien & Fleming boundaries)
# testing median survival 30 vs 20 months in control and treatment group, respectively;
# accrual time 12 and follow-up time 6;
# two-sided significance alpha = 0.05, power 1 - beta = 80%:
lambda1 <- log(2) / 20
lambda2 <- log(2) / 30
allocationRatioPlanned <- 1
accrualTime <- c(0, 12)
followUpTime <- 6
kMax <- 4
alpha <- 0.05
beta <- 0.2
sided <- 2
typeOfDesign <- "OF"

# get the optimal information rates in the group sequential design
informationRates <- getOptimalInformationRates(
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign)
print(informationRates)

# get the group sequential design with the optimal information rates
design <- getDesignGroupSequential(
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  informationRates = informationRates,
  typeOfDesign = typeOfDesign)
summary(design)
print(design)
plot(design)
getDesignCharacteristics(design)

# get the sample sizes in the group sequential design with the optimal information rates
designPlan <- getSampleSizeSurvival(
  design = design,
  lambda1 = lambda1,
  lambda2 = lambda2,
  allocationRatioPlanned = allocationRatioPlanned,
  accrualTime = accrualTime,
  followUpTime = followUpTime)
summary(designPlan)
print(designPlan)
plot(designPlan)

############################################################

# Calculate optimal information rates in a three-stage group sequential design
# (with an equal allocation and O'Brien & Fleming type alpha- and beta-spending boundaries)
# testing H0: pi1 - pi2 = 0 for an alternative H1: pi1 - pi2 = 0.2 with assumed pi2 = 0.2 at event time 12;
# accrual time 12 and follow-up time 6;
# one-sided significance alpha = 0.025, power 1 - beta = 90%:
# thetaH0 <- 1
pi1 <- 0.4
pi2 <- 0.2
allocationRatioPlanned <- 1
eventTime <- 12
accrualTime <- c(0, 12)
followUpTime <- 6
kMax <- 3
alpha <- 0.025
beta <- 0.1
sided <- 1
typeOfDesign <- "asOF"
typeBetaSpending <- "bsOF"

# get the optimal information rates in the group sequential design
informationRates <- getOptimalInformationRates(
  allocationRatioPlanned = allocationRatioPlanned,
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending)
print(informationRates)

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
designPlan <- getSampleSizeSurvival(
  design = design,
  pi1 = pi1,
  pi2 = pi2,
  allocationRatioPlanned = allocationRatioPlanned,
  eventTime = eventTime,
  accrualTime = accrualTime,
  followUpTime = followUpTime)
summary(designPlan)
print(designPlan)
plot(designPlan)

################################################################################
