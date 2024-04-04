#' @title Optimal timing for interim analyses in group sequential trials
#' @author Zhangyi He, Suzie Cro, Laurent Billot

#' OptimInterim (R functions)

################################################################################

#' Get the optimal information rates for interim and final analyses in the group sequential design for testing means or rates in one or two samples
#' Parameter setting
#' @param ... all arguments (starting from "...") can be found in the getDesignGroupSequential and getSampleSizeMeans (or getSampleSizeRates) functions from the rpact package

getOptimalInformationRates <- function(
    ...,
    groups = 2L,
    allocationRatioPlanned = NA_real_,
    kMax = NA_integer_,
    alpha = NA_real_,
    beta = NA_real_,
    sided = 1L,
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
    tolerance = 1e-08,
    seed = NA_real_) {
  options(warn = -1)

  # check and install the optimx, rpact and gtools packages
  if (!require("optimx")) {
    install.packages("optimx")
    library("optimx")
  }
  if (!require("rpact")) {
    install.packages("rpact")
    library("rpact")
  }
  if (!require("gtools")) {
    install.packages("gtools")
    library("gtools")
  }

  # set the seed
  set.seed(ifelse(is.na(seed), 21, seed))

  # define the objective function (the expected sample size under H1) using getSampleSizeMeans
  objectiveFunction <- function(x) {
    informationRates <- sort(unique(c(inv.logit(x), 1))) # 
    if (any(diff(informationRates) < 0.02)) {
      informationRates <- informationRates[-which(diff(informationRates) < 0.02)] # 
    }

    is_NA <- FALSE

    tryCatch(
      {
        design <- getDesignGroupSequential(
          kMax = length(informationRates),
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
        designPlan <- getSampleSizeMeans(
          design = design,
          groups = groups,
          normalApproximation = FALSE,
          meanRatio = FALSE,
          # thetaH0 = ifelse(meanRatio, 1, 0),
          thetaH0 = 0,
          alternative = 0.2,
          stDev = 1,
          allocationRatioPlanned = allocationRatioPlanned)
      }, error = function(e) {
        is_NA <<- TRUE})

    if (is_NA) {
      return(NA)
    } else {
      return(designPlan$expectedNumberOfSubjectsH1[1])
    }
  }

  # # define the objective function (the expected sample size under H1) using getSampleSizeRates
  # objectiveFunction <- function(x) {
  #   informationRates <- sort(unique(c(inv.logit(x), 1))) #
  #   if (any(diff(informationRates) < 0.02)) {
  #     informationRates <- informationRates[-which(diff(informationRates) < 0.02)] #
  #   }
  # 
  #   is_NA <- FALSE
  # 
  #   tryCatch(
  #     {
  #       design <- getDesignGroupSequential(
  #         kMax = length(informationRates),
  #         alpha = alpha,
  #         beta = beta,
  #         sided = sided,
  #         informationRates = informationRates,
  #         futilityBounds = futilityBounds,
  #         typeOfDesign = typeOfDesign,
  #         deltaWT = deltaWT,
  #         deltaPT1 = deltaPT1,
  #         deltaPT0 = deltaPT0,
  #         optimizationCriterion = optimizationCriterion,
  #         gammaA = gammaA,
  #         typeBetaSpending = typeBetaSpending,
  #         userAlphaSpending = userAlphaSpending,
  #         userBetaSpending = userBetaSpending,
  #         gammaB = gammaB,
  #         bindingFutility = bindingFutility,
  #         betaAdjustment = betaAdjustment,
  #         constantBoundsHP = constantBoundsHP,
  #         twoSidedPower = twoSidedPower,
  #         delayedInformation = delayedInformation,
  #         tolerance = tolerance)
  #       designPlan <- getSampleSizeRates(
  #         design = design,
  #         groups = groups,
  #         normalApproximation = TRUE,
  #         riskRatio = FALSE,
  #         # thetaH0 = ifelse(riskRatio, 1, 0),
  #         thetaH0 = 0,
  #         pi1 = 0.4,
  #         pi2 = 0.2,
  #         allocationRatioPlanned = allocationRatioPlanned)
  #     }, error = function(e) {
  #       is_NA <<- TRUE})
  # 
  #   if (is_NA) {
  #     return(NA)
  #   } else {
  #     return(designPlan$expectedNumberOfSubjectsH1[1])
  #   }
  # }
  
  # 
  interimTiming_min <- 1 / kMax
  interimTiming_max <- 1 - 1 / kMax
  interimTiming <- tail(head(seq(0, 1, length.out = kMax + 1), n = -1), n = -1)
  par <- logit(interimTiming)
  res <- optimx(par = par, fn = objectiveFunction, method = c("Nelder-Mead"),
                control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
  # print(res)
  interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
  # print(interimTiming)

  # 
  itn <- 1
  while (any(diff(c(0, interimTiming, 1)) < 0.02) ||
         any(interimTiming < min(interimTiming_min, 0.25)) ||
         any(interimTiming > max(interimTiming_max, 0.95))) {
    itn <- itn + 1

    interimTiming <- seq(max(min(interimTiming), interimTiming_min), 
                         min(max(interimTiming), interimTiming_max), 
                         length.out = kMax - 1) +
      rnorm(n = kMax - 1, mean = 0, sd = 0.05)
    while (any(interimTiming < min(interimTiming_min, 0.25)) || 
           any(interimTiming > max(interimTiming_max, 0.95))) {
      interimTiming <- seq(max(min(interimTiming), interimTiming_min), 
                           min(max(interimTiming), interimTiming_max), 
                           length.out = kMax - 1) +
        rnorm(n = kMax - 1, mean = 0, sd = 0.05)
    }
    par <- logit(interimTiming)
    res <- optimx(par = par, fn = objectiveFunction, method = c("Nelder-Mead"),
                  control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
    # print(res)
    interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
    # print(interimTiming)

    if (itn < 10) {
      next
    } else {
      break
    }
  }

  # 
  itn <- 1
  obj <- ceiling(res$value)
  while (res$value < obj) {
    itn <- itn + 1

    if (any(diff(c(0, interimTiming, 1)) < 0.02) ||
        any(interimTiming < min(interimTiming_min, 0.25)) ||
        any(interimTiming > max(interimTiming_max, 0.95))) {
      interimTiming <- seq(max(min(interimTiming), interimTiming_min), 
                           min(max(interimTiming), interimTiming_max), 
                           length.out = kMax - 1) +
        rnorm(n = kMax - 1, mean = 0, sd = 0.05)
      while (any(interimTiming < min(interimTiming_min, 0.25)) || 
             any(interimTiming > max(interimTiming_max, 0.95))) {
        interimTiming <- seq(max(min(interimTiming), interimTiming_min), 
                             min(max(interimTiming), interimTiming_max), 
                             length.out = kMax - 1) +
          rnorm(n = kMax - 1, mean = 0, sd = 0.05)
      }
      par <- logit(interimTiming)
      res <- optimx(par = par, fn = objectiveFunction, method = c("Nelder-Mead"),
                    control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
      # print(res)
      interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
      # print(interimTiming)
    } else {
      obj <- res$value
      informationRates <- sort(c(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])), 1))

      interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)]))) +
        rnorm(n = kMax - 1, mean = 0, sd = 0.05)
      while (any(interimTiming < min(interimTiming_min, 0.25)) || 
             any(interimTiming > max(interimTiming_max, 0.95))) {
        interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)]))) +
          rnorm(n = kMax - 1, mean = 0, sd = 0.05)
      }
      par <- logit(interimTiming)
      res <- optimx(par = par, fn = objectiveFunction, method = c("Nelder-Mead"),
                    control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
      # print(res)
      interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
      # print(interimTiming)
    }

    if (itn < 10) {
      next
    } else {
      break
    }
  }

  return(informationRates)
}

################################################################################
