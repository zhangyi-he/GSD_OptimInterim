#' @title Optimal timing for interim analyses in group sequential clinical trials
#' @author Zhangyi He, Laurent Billot, Suzie Cro

#' version 1.0

#' rpact plugin (R functions)

################################################################################

#' Get the optimal timing for interim analyses in the group sequential design for continuous endpoints
#' Parameter setting
#' @param weights
#' @param ... all arguments (starting from "...") can be found in the getDesignGroupSequential and getSampleSizeMeans functions from the rpact package

getOptimalInformationRates_ContEndpoint <- function(
    weights = c(0, 0, 1, 0),
    ...,
    groups = 2,
    normalApproximation = FALSE,
    meanRatio = FALSE,
    thetaH0 = ifelse(meanRatio, 1, 0),
    alternative = 0.5,
    stDev = 1,
    allocationRatioPlanned = NA_real_,
    kMax = 3,
    alpha = 0.05,
    beta = 0.1,
    sided = 1,
    futilityBounds = NA_real_,
    typeOfDesign = c("OF", "P", "WT", "PT", "HP", "WToptimum", "asP", "asOF", "asKD", "asHSD", "asUser", "noEarlyEfficacy"),
    deltaWT = NA_real_,
    deltaPT1 = NA_real_,
    deltaPT0 = NA_real_,
    optimizationCriterion = c("ASNH1", "ASNIFH1", "ASNsum"),
    gammaA = NA_real_,
    typeBetaSpending = c("none", "bsP", "bsOF", "bsKD", "bsHSD", "bsUser"),
    userAlphaSpending = NA_real_,
    userBetaSpending = NA_real_,
    gammaB = NA_real_,
    bindingFutility = NA,
    betaAdjustment = NA,
    constantBoundsHP = 3,
    twoSidedPower = NA,
    delayedInformation = NA_real_,
    tolerance = 1e-08) {
  options(warn = -1)

  # check and install the rpact and gtools packages
  if (!require("rpact")) {
    install.packages("rpact")
    library("rpact")
  }
  if (!require("gtools")) {
    install.packages("gtools")
    library("gtools")
  }

  if (is.na(allocationRatioPlanned)) {
    # define the objective function (a weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size)
    fn <- function(x) {
      # convert the input to the information rate
      informationRates <- sort(unique(c(inv.logit(head(x, n = -1)), 1)))
      # if (any(diff(informationRates) < 1e-02)) {
      #   informationRates <- informationRates[-which(diff(informationRates) < 1e-02)]
      # }
      allocationRatioPlanned <- exp(tail(x, n = 1))

      # define the group sequential design
      design <- getDesignGroupSequential(
        kMax = kMax,
        alpha = alpha,
        beta = beta,
        sided = sided,
        informationRates = informationRates,
        futilityBounds = futilityBounds,
        typeOfDesign = typeOfDesign,
        deltaWT = deltaWT,
        deltaPT1 = deltaPT1,
        deltaPT0 = deltaPT0,
        optimizationCriterion = optimizationCriterion,
        gammaA = gammaA,
        typeBetaSpending = typeBetaSpending,
        userAlphaSpending = userAlphaSpending,
        userBetaSpending = userBetaSpending,
        gammaB = gammaB,
        bindingFutility = bindingFutility,
        betaAdjustment = betaAdjustment,
        constantBoundsHP = constantBoundsHP,
        twoSidedPower = twoSidedPower,
        delayedInformation = delayedInformation,
        tolerance = tolerance)

      # produce the sample size (the expected sample sizes under H0, H01, H1 and the maximum sample size)
      designPlan <- getSampleSizeMeans(
        design = design,
        groups = groups,
        normalApproximation = normalApproximation,
        meanRatio = meanRatio,
        thetaH0 = thetaH0,
        alternative = alternative,
        stDev = stDev,
        allocationRatioPlanned = allocationRatioPlanned)
      designPlan <- as.data.frame(designPlan)
      sampleSizes <- c(designPlan$expectedNumberOfSubjectsH0[1],
                       designPlan$expectedNumberOfSubjectsH01[1],
                       designPlan$expectedNumberOfSubjectsH1[1],
                       designPlan$maxNumberOfSubjects[1])

      return(sum(sampleSizes * weights))
    }

    # produce the optimal information rate and the optimal allocation ratio in terms of the weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size
    x <- c(head(tail(seq(0, 1, length.out = kMax + 1), n = -1), n = -1), 1)
    res <- optim(par = c(logit(head(x, n = -1)), log(tail(x, n = 1))), fn = fn, method = "Nelder-Mead")
    informationRates <- c(inv.logit(head(res$par, n = -1)), 1)
    allocationRatioPlanned <- exp(tail(res$par, n = 1))

    return(list("informationRates" = informationRates,
                "allocationRatioPlanned" = allocationRatioPlanned))
  } else {
    # define the objective function (a weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size)
    fn <- function(x) {
      # convert the input to the information rate
      informationRates <- sort(unique(c(inv.logit(x), 1)))
      # if (any(diff(informationRates) < 1e-02)) {
      #   informationRates <- informationRates[-which(diff(informationRates) < 1e-02)]
      # }

      # define the group sequential design
      design <- getDesignGroupSequential(
        kMax = kMax,
        alpha = alpha,
        beta = beta,
        sided = sided,
        informationRates = informationRates,
        futilityBounds = futilityBounds,
        typeOfDesign = typeOfDesign,
        deltaWT = deltaWT,
        deltaPT1 = deltaPT1,
        deltaPT0 = deltaPT0,
        optimizationCriterion = optimizationCriterion,
        gammaA = gammaA,
        typeBetaSpending = typeBetaSpending,
        userAlphaSpending = userAlphaSpending,
        userBetaSpending = userBetaSpending,
        gammaB = gammaB,
        bindingFutility = bindingFutility,
        betaAdjustment = betaAdjustment,
        constantBoundsHP = constantBoundsHP,
        twoSidedPower = twoSidedPower,
        delayedInformation = delayedInformation,
        tolerance = tolerance)

      # produce the sample size (the expected sample sizes under H0, H01, H1 and the maximum sample size)
      designPlan <- getSampleSizeMeans(
        design = design,
        groups = groups,
        normalApproximation = normalApproximation,
        meanRatio = meanRatio,
        thetaH0 = thetaH0,
        alternative = alternative,
        stDev = stDev,
        allocationRatioPlanned = allocationRatioPlanned)
      designPlan <- as.data.frame(designPlan)
      sampleSizes <- c(designPlan$expectedNumberOfSubjectsH0[1],
                       designPlan$expectedNumberOfSubjectsH01[1],
                       designPlan$expectedNumberOfSubjectsH1[1],
                       designPlan$maxNumberOfSubjects[1])

      return(sum(sampleSizes * weights))
    }

    # produce the optimal information rate in terms of the weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size
    x <- head(tail(seq(0, 1, length.out = kMax + 1), n = -1), n = -1)
    res <- optim(par = logit(x), fn = fn, method = "Nelder-Mead")
    informationRates <- c(inv.logit(res$par), 1)

    return(list("informationRates" = informationRates))
  }
}

############################################################

#' Get the optimal timing for interim analyses in group sequential design for binary endpoints
#' Parameter setting
#' @param weights
#' @param ... all arguments (starting from "...") can be found in the getDesignGroupSequential and getSampleSizeRates functions from the rpact package

getOptimalInformationRates_BinEndpoint <- function(
    weights = c(0, 0, 1, 0),
    ...,
    groups = 2,
    normalApproximation = TRUE,
    riskRatio = FALSE,
    thetaH0 = ifelse(riskRatio, 1, 0),
    pi1 = 0.4,
    pi2 = 0.2,
    allocationRatioPlanned = NA_real_,
    kMax = 3,
    alpha = 0.05,
    beta = 0.1,
    sided = 1,
    futilityBounds = NA_real_,
    typeOfDesign = c("OF", "P", "WT", "PT", "HP", "WToptimum", "asP", "asOF", "asKD", "asHSD", "asUser", "noEarlyEfficacy"),
    deltaWT = NA_real_,
    deltaPT1 = NA_real_,
    deltaPT0 = NA_real_,
    optimizationCriterion = c("ASNH1", "ASNIFH1", "ASNsum"),
    gammaA = NA_real_,
    typeBetaSpending = c("none", "bsP", "bsOF", "bsKD", "bsHSD", "bsUser"),
    userAlphaSpending = NA_real_,
    userBetaSpending = NA_real_,
    gammaB = NA_real_,
    bindingFutility = NA,
    betaAdjustment = NA,
    constantBoundsHP = 3,
    twoSidedPower = NA,
    delayedInformation = NA_real_,
    tolerance = 1e-08) {
  options(warn = -1)

  # check and install the rpact and gtools packages
  if (!require("rpact")) {
    install.packages("rpact")
    library("rpact")
  }
  if (!require("gtools")) {
    install.packages("gtools")
    library("gtools")
  }

  if (is.na(allocationRatioPlanned)) {
    # define the objective function (a weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size)
    fn <- function(x) {
      # convert the input to the information rate
      informationRates <- sort(unique(c(inv.logit(head(x, n = -1)), 1)))
      # if (any(diff(informationRates) < 1e-02)) {
      #   informationRates <- informationRates[-which(diff(informationRates) < 1e-02)]
      # }
      allocationRatioPlanned <- exp(tail(x, n = 1))

      # define the group sequential design
      design <- getDesignGroupSequential(
        kMax = kMax,
        alpha = alpha,
        beta = beta,
        sided = sided,
        informationRates = informationRates,
        futilityBounds = futilityBounds,
        typeOfDesign = typeOfDesign,
        deltaWT = deltaWT,
        deltaPT1 = deltaPT1,
        deltaPT0 = deltaPT0,
        optimizationCriterion = optimizationCriterion,
        gammaA = gammaA,
        typeBetaSpending = typeBetaSpending,
        userAlphaSpending = userAlphaSpending,
        userBetaSpending = userBetaSpending,
        gammaB = gammaB,
        bindingFutility = bindingFutility,
        betaAdjustment = betaAdjustment,
        constantBoundsHP = constantBoundsHP,
        twoSidedPower = twoSidedPower,
        delayedInformation = delayedInformation,
        tolerance = tolerance)

      # produce the sample size (the expected sample sizes under H0, H01, H1 and the maximum sample size)
      designPlan <- getSampleSizeMeans(
        design = design,
        groups = groups,
        normalApproximation = normalApproximation,
        riskRatio = riskRatio,
        thetaH0 = thetaH0,
        pi1 = pi1,
        pi2 = pi2,
        allocationRatioPlanned = allocationRatioPlanned)
      designPlan <- as.data.frame(designPlan)
      sampleSizes <- c(designPlan$expectedNumberOfSubjectsH0[1],
                       designPlan$expectedNumberOfSubjectsH01[1],
                       designPlan$expectedNumberOfSubjectsH1[1],
                       designPlan$maxNumberOfSubjects[1])

      return(sum(sampleSizes * weights))
    }

    # produce the optimal information rate and the optimal allocation ratio in terms of the weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size
    x <- c(head(tail(seq(0, 1, length.out = kMax + 1), n = -1), n = -1), 1)
    res <- optim(par = c(logit(head(x, n = -1)), log(tail(x, n = 1))), fn = fn, method = "Nelder-Mead")
    informationRates <- c(inv.logit(head(res$par, n = -1)), 1)
    allocationRatioPlanned <- exp(tail(res$par, n = 1))

    return(list("informationRates" = informationRates,
                "allocationRatioPlanned" = allocationRatioPlanned))
  } else {
    # define the objective function (a weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size)
    fn <- function(x) {
      # convert the input to the information rate
      informationRates <- sort(unique(c(inv.logit(x), 1)))
      # if (any(diff(informationRates) < 1e-02)) {
      #   informationRates <- informationRates[-which(diff(informationRates) < 1e-02)]
      # }

      # define the group sequential design
      design <- getDesignGroupSequential(
        kMax = kMax,
        alpha = alpha,
        beta = beta,
        sided = sided,
        informationRates = informationRates,
        futilityBounds = futilityBounds,
        typeOfDesign = typeOfDesign,
        deltaWT = deltaWT,
        deltaPT1 = deltaPT1,
        deltaPT0 = deltaPT0,
        optimizationCriterion = optimizationCriterion,
        gammaA = gammaA,
        typeBetaSpending = typeBetaSpending,
        userAlphaSpending = userAlphaSpending,
        userBetaSpending = userBetaSpending,
        gammaB = gammaB,
        bindingFutility = bindingFutility,
        betaAdjustment = betaAdjustment,
        constantBoundsHP = constantBoundsHP,
        twoSidedPower = twoSidedPower,
        delayedInformation = delayedInformation,
        tolerance = tolerance)

      # produce the sample size (the expected sample sizes under H0, H01, H1 and the maximum sample size)
      designPlan <- getSampleSizeMeans(
        design = design,
        groups = groups,
        normalApproximation = normalApproximation,
        riskRatio = riskRatio,
        thetaH0 = thetaH0,
        pi1 = pi1,
        pi2 = pi2,
        allocationRatioPlanned = allocationRatioPlanned)
      designPlan <- as.data.frame(designPlan)
      sampleSizes <- c(designPlan$expectedNumberOfSubjectsH0[1],
                       designPlan$expectedNumberOfSubjectsH01[1],
                       designPlan$expectedNumberOfSubjectsH1[1],
                       designPlan$maxNumberOfSubjects[1])

      return(sum(sampleSizes * weights))
    }

    # produce the optimal information rate in terms of the weighted sum of the expected sample sizes under H0, H01, H1 and the maximum sample size
    x <- head(tail(seq(0, 1, length.out = kMax + 1), n = -1), n = -1)
    res <- optim(par = logit(x), fn = fn, method = "Nelder-Mead")
    informationRates <- c(inv.logit(res$par), 1)

    return(list("informationRates" = informationRates))
  }
}

############################################################

#' Get the optimal timing for interim analyses in group sequential design for time-to-event endpoints



################################################################################
