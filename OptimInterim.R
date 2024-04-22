#' OptimInterim

#' R functions

################################################################################

#' Get the optimal information rates for interim and final analyses in the group sequential design for testing means, rates or hazard ratio in two treatment groups.
#' Parameter setting
#' @param ... all arguments (starting from "...") can be found in the getDesignGroupSequential and getSampleSizeMeans (or getSampleSizeRates or getSampleSizeSurvival) functions from the rpact package

getOptimalInformationRates <- function(
    ...,
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

  if (kMax > 1) { # a group sequential design
    # calculate the objective (the expected sample size under H1)
    # using getDesignGroupSequential and getSampleSizeMeans (or getSampleSizeRates or getSampleSizeSurvival)
    calculateObjective <- function(x) {
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
            groups = 2,
            normalApproximation = FALSE,
            meanRatio = FALSE,
            # thetaH0 = ifelse(meanRatio, 1, 0),
            thetaH0 = 0,
            alternative = 0.2,
            stDev = 1,
            allocationRatioPlanned = allocationRatioPlanned)
          # designPlan <- getSampleSizeRates(
          #   design = design,
          #   groups = 2,
          #   normalApproximation = TRUE,
          #   riskRatio = FALSE,
          #   # thetaH0 = ifelse(riskRatio, 1, 0),
          #   thetaH0 = 0,
          #   pi1 = 0.4,
          #   pi2 = 0.2,
          #   allocationRatioPlanned = allocationRatioPlanned)
          # designPlan <- getSampleSizeSurvival(
          #   design = design,
          #   thetaH0 = 1,
          #   pi1 = 0.2,
          #   pi2 = 0.4,
          #   allocationRatioPlanned = allocationRatioPlanned,
          #   eventTime = 12,
          #   accrualTime = c(0, 12),
          #   followUpTime = 6)
        }, error = function(e) {
          is_NA <<- TRUE})

      if (is_NA) {
        return(NA)
      } else {
        return(designPlan$expectedNumberOfSubjectsH1[1])
        # return(designPlan$expectedEventsH1[1])
      }
    }

    # calculate the rejection (the early rejection probability per stage)
    # using getDesignGroupSequential and getSampleSizeMeans (or getSampleSizeRates or getSampleSizeSurvival)
    calculateRejection <- function(t) {
      informationRates <- sort(c(t, 1)) #

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
            groups = 2,
            normalApproximation = FALSE,
            meanRatio = FALSE,
            # thetaH0 = ifelse(meanRatio, 1, 0),
            thetaH0 = 0,
            alternative = 0.2,
            stDev = 1,
            allocationRatioPlanned = allocationRatioPlanned)
          # designPlan <- getSampleSizeRates(
          #   design = design,
          #   groups = 2,
          #   normalApproximation = TRUE,
          #   riskRatio = FALSE,
          #   # thetaH0 = ifelse(riskRatio, 1, 0),
          #   thetaH0 = 0,
          #   pi1 = 0.4,
          #   pi2 = 0.2,
          #   allocationRatioPlanned = allocationRatioPlanned)
          # designPlan <- getSampleSizeSurvival(
          #   design = design,
          #   thetaH0 = 1,
          #   pi1 = 0.2,
          #   pi2 = 0.4,
          #   allocationRatioPlanned = allocationRatioPlanned,
          #   eventTime = 12,
          #   accrualTime = c(0, 12),
          #   followUpTime = 6)
        }, error = function(e) {
          is_NA <<- TRUE})

      if (is_NA) {
        return(NA)
      } else {
        return(head(as.numeric(designPlan$rejectPerStage), n = -1))
      }
    }

    #
    interimTiming_min <- min(1 / kMax, 0.25)
    # interimTiming_min <- 0.01
    # continue <- TRUE
    # is_NA <- FALSE
    # while (continue) {
    #   tryCatch(
    #     {
    #       design <- getDesignGroupSequential(
    #         kMax = 2,
    #         alpha = alpha,
    #         beta = beta,
    #         sided = sided,
    #         informationRates = sort(c(interimTiming_min, 1)),
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
    #       designPlan <- getSampleSizeMeans(
    #         design = design,
    #         groups = 2,
    #         normalApproximation = FALSE,
    #         meanRatio = FALSE,
    #         # thetaH0 = ifelse(meanRatio, 1, 0),
    #         thetaH0 = 0,
    #         alternative = 0.2,
    #         stDev = 1,
    #         allocationRatioPlanned = allocationRatioPlanned)
    #       # designPlan <- getSampleSizeRates(
    #       #   design = design,
    #       #   groups = 2,
    #       #   normalApproximation = TRUE,
    #       #   riskRatio = FALSE,
    #       #   # thetaH0 = ifelse(riskRatio, 1, 0),
    #       #   thetaH0 = 0,
    #       #   pi1 = 0.4,
    #       #   pi2 = 0.2,
    #       #   allocationRatioPlanned = allocationRatioPlanned)
    #       # designPlan <- getSampleSizeSurvival(
    #       #   design = design,
    #       #   thetaH0 = 1,
    #       #   pi1 = 0.2,
    #       #   pi2 = 0.4,
    #       #   allocationRatioPlanned = allocationRatioPlanned,
    #       #   eventTime = 12,
    #       #   accrualTime = c(0, 12),
    #       #   followUpTime = 6)
    #       designPlan <- as.data.frame(designPlan)
    #     }, error = function(e) {
    #       is_NA <<- TRUE})

    #   if (is_NA) {
    #     continue <- TRUE
    #     interimTiming_min <- interimTiming_min + 0.01
    #     is_NA <- FALSE
    #   } else {
    #     continue <-
    #       (abs(designPlan$expectedNumberOfSubjectsH1[1] - designPlan$maxNumberOfSubjects[1]) < 1)
    #     # continue <-
    #     #   (abs(designPlan$expectedEventsH1[1] - designPlan$maxNumberOfEvents[1]) < 1)
    #     interimTiming_min <- ifelse(continue, interimTiming_min + 0.01, interimTiming_min)
    #   }
    # }
    # # print(interimTiming_min)

    #
    interimTiming_max <- max(1 - 1 / kMax, 0.95)
    # interimTiming_max <- 0.99
    # continue <- TRUE
    # is_NA <- FALSE
    # while (continue) {
    #   tryCatch(
    #     {
    #       design <- getDesignGroupSequential(
    #         kMax = 2,
    #         alpha = alpha,
    #         beta = beta,
    #         sided = sided,
    #         informationRates = sort(c(interimTiming_max, 1)),
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
    #       designPlan <- getSampleSizeMeans(
    #         design = design,
    #         groups = 2,
    #         normalApproximation = FALSE,
    #         meanRatio = FALSE,
    #         # thetaH0 = ifelse(meanRatio, 1, 0),
    #         thetaH0 = 0,
    #         alternative = 0.2,
    #         stDev = 1,
    #         allocationRatioPlanned = allocationRatioPlanned)
    #       # designPlan <- getSampleSizeRates(
    #       #   design = design,
    #       #   groups = 2,
    #       #   normalApproximation = TRUE,
    #       #   riskRatio = FALSE,
    #       #   # thetaH0 = ifelse(riskRatio, 1, 0),
    #       #   thetaH0 = 0,
    #       #   pi1 = 0.4,
    #       #   pi2 = 0.2,
    #       #   allocationRatioPlanned = allocationRatioPlanned)
    #       # designPlan <- getSampleSizeSurvival(
    #       #   design = design,
    #       #   thetaH0 = 1,
    #       #   pi1 = 0.2,
    #       #   pi2 = 0.4,
    #       #   allocationRatioPlanned = allocationRatioPlanned,
    #       #   eventTime = 12,
    #       #   accrualTime = c(0, 12),
    #       #   followUpTime = 6)
    #       designPlan <- as.data.frame(designPlan)
    #     }, error = function(e) {
    #       is_NA <<- TRUE})

    #   if (is_NA) {
    #     continue <- TRUE
    #     interimTiming_max <- interimTiming_max - 0.01
    #     is_NA <- FALSE
    #   } else {
    #     continue <-
    #       (abs(designPlan$expectedNumberOfSubjectsH1[1] - designPlan$maxNumberOfSubjects[1]) < 1)
    #     # continue <-
    #     #   (abs(designPlan$expectedEventsH1[1] - designPlan$maxNumberOfEvents[1]) < 1)
    #     interimTiming_max <- ifelse(continue, interimTiming_max - 0.01, interimTiming_max)
    #   }
    # }
    # print(interimTiming_max)

    # initialise optimisation
    message("Initialising optimisation...")
    itn <- 1
    # print(paste0("iteration ", itn, ": initialise and run optimisation"))

    # initialise the input
    interimTiming <- tail(head(seq(from = 0, to = 1, length.out = kMax + 1), n = -1), n = -1)
    # interimTiming <- tail(head(seq(from = interimTiming_min,
    #                                to = interimTiming_max,
    #                                length.out = kMax + 1), n = -1), n = -1)
    par <- logit(interimTiming)

    # optimise the objective
    res <- optimx(par = par, fn = calculateObjective, method = c("Nelder-Mead"),
                  control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
    # print(res)
    interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
    # print(interimTiming)
    rejectProb <- calculateRejection(interimTiming)
    # print(rejectProb)

    # check and update optimisation
    while (any(diff(c(0, interimTiming, 1)) < 0.02) ||
           any(rejectProb < 0.01)) {
      itn <- itn + 1
      # print(paste0("iteration ", itn, ": check and update optimisation"))

      # initialise the input
      if (any(rejectProb < 0.01)) {
        interimTiming <- interimTiming[-which(rejectProb < 0.01)]
        interimTiming <- seq(from = max(min(interimTiming), interimTiming_min),
                             to = min(max(interimTiming), interimTiming_max),
                             length.out = kMax - 1) +
          rnorm(n = kMax - 1, mean = 0, sd = 0.05)
      } else {
        interimTiming <- seq(from = max(min(interimTiming), interimTiming_min),
                             to = min(max(interimTiming), interimTiming_max),
                             length.out = kMax - 1) +
          rnorm(n = kMax - 1, mean = 0, sd = 0.05)
      }
      interimTiming <- sort(interimTiming)
      interimTiming[1] <- max(min(interimTiming), interimTiming_min)
      interimTiming[kMax - 1] <- min(max(interimTiming), interimTiming_max)
      if (any(diff(c(0, interimTiming, 1)) < 0.02)) {
        par <- logit(seq(from = min(interimTiming),
                         to = max(interimTiming),
                         length.out = kMax - 1))
      } else {
        par <- logit(interimTiming)
      }

      # optimise the objective
      res <- optimx(par = par, fn = calculateObjective, method = c("Nelder-Mead"),
                    control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
      # print(res)
      interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
      # print(interimTiming)
      rejectProb <- calculateRejection(interimTiming)
      # print(rejectProb)

      if (itn < 10) {
        next
      } else {
        break
      }
    }

    # check and update optimisation
    message("Updating optimisation...")
    itn <- 1
    obj <- ceiling(res$value)
    while (res$value < obj) {
      itn <- itn + 1
      # print(paste0("iteration ", itn, ": check and update optimisation"))

      if (any(diff(c(0, interimTiming, 1)) < 0.02) ||
          any(rejectProb < 0.01)) {
        # initialise the input
        if (any(rejectProb < 0.01)) {
          interimTiming <- interimTiming[-which(rejectProb < 0.01)]
          interimTiming <- seq(from = max(min(interimTiming), interimTiming_min),
                               to = min(max(interimTiming), interimTiming_max),
                               length.out = kMax - 1) +
            rnorm(n = kMax - 1, mean = 0, sd = 0.05)
        } else {
          interimTiming <- seq(from = max(min(interimTiming), interimTiming_min),
                               to = min(max(interimTiming), interimTiming_max),
                               length.out = kMax - 1) +
            rnorm(n = kMax - 1, mean = 0, sd = 0.05)
        }
        interimTiming <- sort(interimTiming)
        interimTiming[1] <- max(min(interimTiming), interimTiming_min)
        interimTiming[kMax - 1] <- min(max(interimTiming), interimTiming_max)
        if (any(diff(c(0, interimTiming, 1)) < 0.02)) {
          par <- logit(seq(from = min(interimTiming),
                           to = max(interimTiming),
                           length.out = kMax - 1))
        } else {
          par <- logit(interimTiming)
        }

        # optimise the objective
        res <- optimx(par = par, fn = calculateObjective, method = c("Nelder-Mead"),
                      control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
        # print(res)
        interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
        # print(interimTiming)
        rejectProb <- calculateRejection(interimTiming)
        # print(rejectProb)
      } else {
        obj <- res$value
        informationRates <- sort(c(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])), 1))

        # initialise the input
        interimTiming <- inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])) +
          rnorm(n = kMax - 1, mean = 0, sd = 0.05)
        interimTiming <- sort(interimTiming)
        interimTiming[1] <- max(min(interimTiming), interimTiming_min)
        interimTiming[kMax - 1] <- min(max(interimTiming), interimTiming_max)
        if (any(diff(c(0, interimTiming, 1)) < 0.02)) {
          par <- logit(seq(from = min(interimTiming),
                           to = max(interimTiming),
                           length.out = kMax - 1))
        } else {
          par <- logit(interimTiming)
        }

        # optimise the objective
        res <- optimx(par = par, fn = calculateObjective, method = c("Nelder-Mead"),
                      control = list(trace = 0, follow.on = TRUE, maximize = FALSE, maxit = 1e+04))
        # print(res)
        interimTiming <- sort(inv.logit(as.numeric(res["Nelder-Mead", 1:(kMax - 1)])))
        # print(interimTiming)
        rejectProb <- calculateRejection(interimTiming)
        # print(rejectProb)
      }

      if (itn < 10) {
        next
      } else {
        break
      }
    }

    # check the output
    if (any(diff(c(0, informationRates)) < 0.02) ||
        any(rejectProb < 0.01)) {
      informationRates <- NA
      message("Failure in optimisation!")
    } else {
      message("Success in optimisation!")
    }
  } else { # a fixed sample design
    informationRates <- 1
  }

  return(informationRates)
}

################################################################################
